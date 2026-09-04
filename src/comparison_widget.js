/**
 * OpenFE Benchmark Comparison Widget
 *
 * Interactive d3-based visualization for multi-force-field benchmark comparisons.
 * Renders grouped bar charts, ECDF curves, scatter plots, and 2D structure depictions.
 */

(function() {
  'use strict';

  // ============================================================================
  // Constants and Configuration
  // ============================================================================

  const COLORS = {
    purple: '#8b5cf6',
    blue: '#3b82f6',
    red: '#ef4444',
    gray: '#9ca3af',
  };

  const SERIES_PALETTE = [
    '#8b5cf6',
    '#ef4444',
    '#10b981',
    '#f59e0b',
    '#ec4899',
    '#84cc16',
    '#f97316',
    '#e11d48',
    '#a16207',
    '#7c3aed',
  ];

  const PANEL_METRICS = [
    { id: 'ddg_mue', label: 'ΔΔG MUE', ylabel: 'MUE (kcal/mol)' },
    { id: 'ddg_rmse', label: 'ΔΔG RMSE', ylabel: 'RMSE (kcal/mol)' },
    { id: 'dg_rmse', label: 'ΔG RMSE (centralized)', ylabel: 'RMSE (kcal/mol)' },
    { id: 'dg_ktau', label: 'ΔG Kendall τ', ylabel: 'Kendall τ' }
  ];

  // ============================================================================
  // Global State
  // ============================================================================

  let data = null;
  let colorScale = null;
  let selectedEdge = null; // { edge_id, ligand_a, ligand_b, system_id } or null
  let selectedLigand = null; // { ligand, system_id } or null
  let activeDetailSystemId = null;
  let detailPanelResizeObserver = null;
  let detailPanelLastWidth = null;
  let detailPanelResizeRaf = null;

  // ============================================================================
  // Data Loading and Setup
  // ============================================================================

  function loadData() {
    const dataEl = document.getElementById('widget-data');
    if (!dataEl) {
      throw new Error('Data element #widget-data not found');
    }
    return JSON.parse(dataEl.textContent);
  }

  function setupColorScale(keys) {
    return d3.scaleOrdinal()
      .domain(keys)
      .range(SERIES_PALETTE);
  }

  function renderLatex(text) {
    if (typeof katex !== 'undefined') {
      try {
        // Extract LaTeX expressions and render them
        const regex = /\$([^$]+)\$/g;
        return text.replace(regex, (match, latex) => {
          try {
            return katex.renderToString(latex, { throwOnError: false });
          } catch (e) {
            return match;
          }
        });
      } catch (e) {
        return text;
      }
    }
    // KaTeX not available - return raw text as-is
    return text;
  }

  // ============================================================================
  // Main Panel: Grouped Bar Charts
  // ============================================================================

  function renderMainPanel(container) {
    const width = 1200;
    const height = 800;

    // Calculate dynamic top margin based on legend height
    // Legend: title (15px) + keys (25px each) + buffer (30px)
    const legendHeight = 15 + (data.keys.length * 25) + 30;
    const topMargin = Math.max(140, legendHeight);

    const margin = { top: topMargin, right: 120, bottom: 80, left: 80 };
    const panelWidth = (width - margin.left - margin.right - 80) / 2;
    const panelHeight = (height - margin.top - margin.bottom - 80) / 2;

    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .attr('class', 'main-panel');

    // Render 2x2 grid of bar charts
    PANEL_METRICS.forEach((metric, idx) => {
      const row = Math.floor(idx / 2);
      const col = idx % 2;
      const x = margin.left + col * (panelWidth + 60);
      const y = margin.top + row * (panelHeight + 60);

      renderBarChart(svg, metric, x, y, panelWidth, panelHeight);
    });

    // Shared legend (top left)
    renderLegend(svg, 10, 10);

    svg.append('text')
      .attr('x', width / 2)
      .attr('y', height - 32)
      .attr('text-anchor', 'middle')
      .attr('font-size', '13px')
      .attr('fill', '#4b5563')
      .text('Click a system name on the x-axis to open its detail panel.');

    if (data.reference_model) {
      svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 12)
        .attr('text-anchor', 'middle')
        .attr('font-size', '12px')
        .attr('fill', '#6b7280')
        .text(`Hatched bars differ significantly from ${data.reference_model}.`);
    }
  }

  function renderBarChart(svg, metric, x, y, width, height) {
    const g = svg.append('g')
      .attr('transform', `translate(${x},${y})`)
      .attr('class', `bar-chart bar-chart-${metric.id}`);

    // Define a single hatch pattern for bars that differ from the reference model.
    const defs = svg.append('defs');
    defs.append('pattern')
      .attr('id', `hatch-different-${metric.id}`)
      .attr('patternUnits', 'userSpaceOnUse')
      .attr('width', 6)
      .attr('height', 6)
      .append('path')
      .attr('d', 'M-1,1 l2,-2 M0,6 l6,-6 M5,7 l2,-2')
      .attr('stroke', '#000')
      .attr('stroke-width', 0.75);

    // Title
    g.append('text')
      .attr('x', width / 2)
      .attr('y', -10)
      .attr('text-anchor', 'middle')
      .attr('font-weight', 'bold')
      .attr('font-size', '14px')
      .text(metric.label);

    // Collect data
    const chartData = [];
    data.systems.forEach(sys => {
      if (data.bars[sys.id] && data.bars[sys.id][metric.id]) {
        data.bars[sys.id][metric.id].forEach(bar => {
          chartData.push({
            system: sys.system_name,
            system_id: sys.id,
            key: bar.key,
            value: bar.value,
            ci_lower: bar.ci_lower,
            ci_upper: bar.ci_upper,
            cld: bar.cld,
            different_from_reference: Boolean(bar.different_from_reference)
          });
        });
      }
    });

    if (chartData.length === 0) {
      g.append('text')
        .attr('x', width / 2)
        .attr('y', height / 2)
        .attr('text-anchor', 'middle')
        .attr('fill', '#999')
        .text('No data');
      return;
    }

    // Scales
    const systems = [...new Set(chartData.map(d => d.system))];
    const x0 = d3.scaleBand()
      .domain(systems)
      .range([0, width])
      .padding(0.2);

    const x1 = d3.scaleBand()
      .domain(data.keys)
      .range([0, x0.bandwidth()])
      .padding(0.05);

    const maxVal = d3.max(chartData, d => d.ci_upper) || 1;
    const minVal = metric.id === 'dg_ktau'
      ? Math.min(0, d3.min(chartData, d => d.ci_lower) || 0)
      : 0;

    const yScale = d3.scaleLinear()
      .domain([minVal, maxVal * 1.1])
      .range([height, 0])
      .nice();

    // Axes
    const xAxis = d3.axisBottom(x0);
    const yAxis = d3.axisLeft(yScale).ticks(5);

    g.append('g')
      .attr('class', 'x-axis')
      .attr('transform', `translate(0,${height})`)
      .call(xAxis)
      .selectAll('text')
      .attr('class', 'x-axis-label')
      .attr('cursor', 'pointer')
      .attr('tabindex', 0)
      .style('font-weight', 'bold')
      .on('click', function(event, d) {
        const systemId = data.systems.find(s => s.system_name === d)?.id;
        if (systemId) openDetailPanel(systemId);
      })
      .on('mouseenter', function() {
        d3.select(this).style('fill', '#3b82f6');
      })
      .on('mouseleave', function() {
        d3.select(this).style('fill', 'black');
      })
      .on('keydown', function(event, d) {
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          const systemId = data.systems.find(s => s.system_name === d)?.id;
          if (systemId) openDetailPanel(systemId);
        }
      });

    g.append('g')
      .attr('class', 'y-axis')
      .call(yAxis);

    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -height / 2)
      .attr('y', -50)
      .attr('text-anchor', 'middle')
      .attr('font-size', '12px')
      .text(metric.ylabel);

    // Bars
    const barGroups = g.selectAll('.bar-group')
      .data(systems)
      .join('g')
      .attr('class', 'bar-group')
      .attr('transform', d => `translate(${x0(d)},0)`);

    barGroups.each(function(system) {
      const group = d3.select(this);
      const systemData = chartData.filter(d => d.system === system);

      systemData.forEach(bar => {
        const barG = group.append('g')
          .attr('class', `bar-item bar-${bar.key.replace(/\s+/g, '-')}`);

        // Bar rect with solid color fill
        const barHeight = Math.max(0, height - yScale(bar.value));
        const barRect = barG.append('rect')
          .attr('x', x1(bar.key))
          .attr('y', yScale(bar.value))
          .attr('width', x1.bandwidth())
          .attr('height', barHeight)
          .attr('fill', colorScale(bar.key))
          .attr('stroke', 'black')
          .attr('stroke-width', 0.5);

        if (bar.different_from_reference) {
          barG.append('rect')
            .attr('x', x1(bar.key))
            .attr('y', yScale(bar.value))
            .attr('width', x1.bandwidth())
            .attr('height', barHeight)
            .attr('fill', `url(#hatch-different-${metric.id})`)
            .attr('stroke', 'none')
            .style('pointer-events', 'none');
        }

        // Error bars
        const errorX = x1(bar.key) + x1.bandwidth() / 2;
        barG.append('line')
          .attr('x1', errorX)
          .attr('x2', errorX)
          .attr('y1', yScale(bar.ci_lower))
          .attr('y2', yScale(bar.ci_upper))
          .attr('stroke', 'black')
          .attr('stroke-width', 1.5);

        barG.append('line')
          .attr('x1', errorX - 3)
          .attr('x2', errorX + 3)
          .attr('y1', yScale(bar.ci_lower))
          .attr('y2', yScale(bar.ci_lower))
          .attr('stroke', 'black')
          .attr('stroke-width', 1.5);

        barG.append('line')
          .attr('x1', errorX - 3)
          .attr('x2', errorX + 3)
          .attr('y1', yScale(bar.ci_upper))
          .attr('y2', yScale(bar.ci_upper))
          .attr('stroke', 'black')
          .attr('stroke-width', 1.5);
      });
    });
  }

  function renderLegend(svg, x, y) {
    const legend = svg.append('g')
      .attr('class', 'legend')
      .attr('transform', `translate(${x},${y})`);

    legend.append('text')
      .attr('x', 0)
      .attr('y', 0)
      .attr('font-weight', 'bold')
      .attr('font-size', '12px')
      .text('Force Fields:');

    data.keys.forEach((key, i) => {
      const item = legend.append('g')
        .attr('transform', `translate(0,${20 + i * 25})`);

      item.append('rect')
        .attr('width', 18)
        .attr('height', 18)
        .attr('fill', colorScale(key));

      item.append('text')
        .attr('x', 25)
        .attr('y', 14)
        .attr('font-size', '11px')
        .attr('fill', 'black')
        .text(key.replace(/\$/g, ''));
    });
  }

  // ============================================================================
  // Detail Panel
  // ============================================================================

  function openDetailPanel(systemId) {
    const system = data.systems.find(s => s.id === systemId);
    if (!system) return;

    // Clear and show detail panel
    const detailPanel = d3.select('#detail-panel');
    if (detailPanel.empty()) {
      const container = d3.select('#widget-root')
        .append('div')
        .attr('id', 'detail-panel')
        .style('position', 'fixed')
        .style('right', '0')
        .style('top', '0')
        .style('width', '50%')
        .style('height', '100vh')
        .style('background', 'white')
        .style('border-left', '2px solid #ccc')
        .style('overflow-y', 'auto')
        .style('padding', '20px')
        .style('box-shadow', '-2px 0 5px rgba(0,0,0,0.1)');

      const closeBtn = container.append('button')
        .attr('class', 'close-btn')
        .style('position', 'absolute')
        .style('top', '10px')
        .style('right', '10px')
        .style('padding', '5px 10px')
        .style('cursor', 'pointer')
        .text('Close')
        .on('click', () => {
          if (detailPanelResizeObserver) {
            detailPanelResizeObserver.disconnect();
            detailPanelResizeObserver = null;
          }
          if (detailPanelResizeRaf !== null) {
            window.cancelAnimationFrame(detailPanelResizeRaf);
            detailPanelResizeRaf = null;
          }
          activeDetailSystemId = null;
          detailPanelLastWidth = null;
          container.remove();
        });

      setupDetailPanelResizeBehavior(container.node());
    }

    activeDetailSystemId = systemId;
    renderDetailContent(systemId);
  }

  function setupDetailPanelResizeBehavior(panelEl) {
    if (typeof ResizeObserver === 'undefined') {
      return;
    }

    if (detailPanelResizeObserver) {
      detailPanelResizeObserver.disconnect();
    }

    detailPanelResizeObserver = new ResizeObserver(entries => {
      const entry = entries[0];
      if (!entry || !activeDetailSystemId) return;

      const width = Math.round(entry.contentRect.width);
      if (width === detailPanelLastWidth) {
        return;
      }
      detailPanelLastWidth = width;

      if (detailPanelResizeRaf !== null) {
        window.cancelAnimationFrame(detailPanelResizeRaf);
      }
      detailPanelResizeRaf = window.requestAnimationFrame(() => {
        detailPanelResizeRaf = null;
        if (activeDetailSystemId && document.getElementById('detail-panel')) {
          renderDetailContent(activeDetailSystemId);
        }
      });
    });

    detailPanelResizeObserver.observe(panelEl);
    detailPanelLastWidth = Math.round(panelEl.getBoundingClientRect().width);
  }

  function renderDetailContent(systemId) {
    const container = d3.select('#detail-panel');
    container.selectAll('*:not(.close-btn)').remove();

    const system = data.systems.find(s => s.id === systemId);

    // Title
    container.append('h2')
      .text(system.system_name);

    // Missing/dropped notes
    if (data.missing[systemId] && data.missing[systemId].length > 0) {
      container.append('p')
        .style('color', '#f59e0b')
        .html(`<strong>No data for:</strong> ${data.missing[systemId].join(', ')}`);
    }

    if (data.dropped[systemId]) {
      const dropped = data.dropped[systemId];
      const edgeCount = (dropped.edges || []).length;
      const nodeCount = (dropped.nodes || []).length;
      if (edgeCount > 0 || nodeCount > 0) {
        const parts = [];
        if (edgeCount > 0) parts.push(`${edgeCount} edge${edgeCount !== 1 ? 's' : ''}`);
        if (nodeCount > 0) parts.push(`${nodeCount} node${nodeCount !== 1 ? 's' : ''}`);
        container.append('p')
          .style('color', '#ef4444')
          .html(`<strong>Excluded from cross-key statistics:</strong> ${parts.join(', ')}`);
      }
    }

    // Structure strip (initially empty)
    const structureStrip = container.append('div')
      .attr('id', 'structure-strip')
      .style('margin', '20px 0')
      .style('padding', '15px')
      .style('border', '1px solid #e5e7eb')
      .style('border-radius', '4px')
      .style('background', '#f9fafb');

    structureStrip.append('div')
      .attr('class', 'structure-placeholder')
      .style('color', '#9ca3af')
      .style('font-style', 'italic')
      .text('Click a scatter point to view 2D structures');

    // ECDF
    container.append('h3').text('ECDF (Cumulative Error Distribution)');
    const ecdfDiv = container.append('div').attr('id', 'ecdf-chart');
    renderECDF(ecdfDiv.node(), systemId);

    // Scatters
    container.append('h3').text('Per-Force-Field Scatters (click one of the points)');
    const scatterDiv = container.append('div').attr('id', 'scatter-grid');
    renderScatterGrid(scatterDiv.node(), systemId);

    // Pair plot control
    container.append('h3').text('Key-vs-Key Pair Plot');
    const pairDiv = container.append('div')
      .attr('id', 'pair-plot-control')
      .style('margin-bottom', '40px');
    renderPairPlotControl(pairDiv.node(), systemId);
  }

  // ============================================================================
  // ECDF with Hover Readout
  // ============================================================================

  function renderECDF(container, systemId) {
    if (!data.ecdf[systemId]) {
      d3.select(container).append('p').text('No ECDF data');
      return;
    }

    const availableWidth = getContainerWidth(container, 600);
    const width = Math.max(520, Math.min(980, availableWidth));
    const height = Math.max(340, Math.min(560, Math.round(width * 0.62)));
    const margin = { top: 20, right: 120, bottom: 50, left: 60 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height);

    const g = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // Scales
    const ecdfData = data.ecdf[systemId];
    const allX = [];
    Object.values(ecdfData).forEach(series => {
      allX.push(...series.x);
    });

    const xScale = d3.scaleLinear()
      .domain([0, d3.max(allX) || 1])
      .range([0, plotWidth])
      .nice();

    const yScale = d3.scaleLinear()
      .domain([0, 1])
      .range([plotHeight, 0]);

    // Axes
    g.append('g')
      .attr('transform', `translate(0,${plotHeight})`)
      .call(d3.axisBottom(xScale));

    g.append('g')
      .call(d3.axisLeft(yScale).ticks(5));

    g.append('text')
      .attr('x', plotWidth / 2)
      .attr('y', height - 32)
      .attr('text-anchor', 'middle')
      .text('|ΔΔG calc - ΔΔG exp| (kcal/mol)');

    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -plotHeight / 2)
      .attr('y', -45)
      .attr('text-anchor', 'middle')
      .text('Cumulative Probability');

    // Draw curves
    data.keys.forEach(key => {
      if (!ecdfData[key]) return;

      const series = ecdfData[key];
      const line = d3.line()
        .x((d, i) => xScale(series.x[i]))
        .y((d, i) => yScale(series.y[i]))
        .curve(d3.curveStepAfter);

      g.append('path')
        .datum(series.y)
        .attr('class', `ecdf-curve ecdf-${key.replace(/\s+/g, '-')}`)
        .attr('d', line)
        .attr('fill', 'none')
        .attr('stroke', colorScale(key))
        .attr('stroke-width', 2);

      // CI bands if available
      if (series.ci_lower && series.ci_upper) {
        const area = d3.area()
          .x((d, i) => xScale(series.x[i]))
          .y0((d, i) => yScale(series.ci_lower[i]))
          .y1((d, i) => yScale(series.ci_upper[i]))
          .curve(d3.curveStepAfter);

        g.append('path')
          .datum(series.y)
          .attr('class', `ecdf-ci ecdf-ci-${key.replace(/\s+/g, '-')}`)
          .attr('d', area)
          .attr('fill', colorScale(key))
          .attr('fill-opacity', 0.2)
          .attr('stroke', 'none');
      }
    });

    // Hover readout
    setupECDFHover(g, ecdfData, xScale, yScale, plotWidth, plotHeight);
  }

  function setupECDFHover(g, ecdfData, xScale, yScale, width, height) {
    const hoverG = g.append('g')
      .attr('class', 'ecdf-hover')
      .style('display', 'none');

    const verticalLine = hoverG.append('line')
      .attr('class', 'hover-vertical')
      .attr('y1', 0)
      .attr('y2', height)
      .attr('stroke', '#999')
      .attr('stroke-width', 1)
      .attr('stroke-dasharray', '3,3');

    const xValueLabel = hoverG.append('text')
      .attr('class', 'hover-x-label')
      .attr('y', height + 15)
      .attr('text-anchor', 'middle')
      .attr('font-size', '11px')
      .attr('fill', '#333');

    const horizontalLinesG = hoverG.append('g').attr('class', 'hover-horizontal-lines');

    const overlay = g.append('rect')
      .attr('class', 'ecdf-overlay')
      .attr('width', width)
      .attr('height', height)
      .attr('fill', 'none')
      .attr('pointer-events', 'all');

    overlay.on('mousemove', function(event) {
      const [mouseX] = d3.pointer(event);
      const xValue = xScale.invert(mouseX);

      hoverG.style('display', null);
      verticalLine.attr('x1', mouseX).attr('x2', mouseX);
      xValueLabel.attr('x', mouseX).text(xValue.toFixed(2));

      // Compute y values for each series using step semantics
      const readings = [];
      data.keys.forEach(key => {
        if (!ecdfData[key]) return;

        const series = ecdfData[key];
        const idx = d3.bisectRight(series.x, xValue);
        const yValue = idx === 0 ? 0 : (idx >= series.x.length ? 1 : series.y[idx - 1]);

        readings.push({
          key,
          yValue,
          yPixel: yScale(yValue),
          color: colorScale(key)
        });
      });

      // De-overlap labels
      const minSpacing = 18;
      readings.sort((a, b) => a.yValue - b.yValue);

      for (let i = 0; i < readings.length; i++) {
        if (i > 0) {
          const prev = readings[i - 1];
          const minY = prev.yPixel - minSpacing;
          if (readings[i].yPixel > minY) {
            readings[i].yPixel = minY;
          }
        }
        readings[i].yPixel = Math.max(0, Math.min(height, readings[i].yPixel));
      }

      // Render horizontal lines and labels
      horizontalLinesG.selectAll('*').remove();

      readings.forEach(r => {
        const lineG = horizontalLinesG.append('g');

        lineG.append('line')
          .attr('x1', 0)
          .attr('x2', xScale(xValue))
          .attr('y1', yScale(r.yValue))
          .attr('y2', yScale(r.yValue))
          .attr('stroke', r.color)
          .attr('stroke-width', 1)
          .attr('stroke-dasharray', '2,2');

        lineG.append('text')
          .attr('x', -5)
          .attr('y', r.yPixel)
          .attr('dy', '0.32em')
          .attr('text-anchor', 'end')
          .attr('font-size', '10px')
          .attr('fill', r.color)
          .text(r.yValue.toFixed(2));
      });
    });

    overlay.on('mouseleave', function() {
      hoverG.style('display', 'none');
    });
  }

  // ============================================================================
  // Scatter Grids
  // ============================================================================

  function renderScatterGrid(container, systemId) {
    const scatterDiv = d3.select(container);
    const layout = computeScatterLayout(container, data.keys.length);

    // Calculate global axis limits across ALL systems for this detail view
    let allDDG = [];
    let allDG = [];

    // Collect all ΔΔG values
    Object.values(data.edges).forEach(systemEdges => {
      systemEdges.forEach(edge => {
        allDDG.push(edge.ddg_exp);
        data.keys.forEach(key => {
          if (edge.per_key[key]) {
            allDDG.push(edge.per_key[key].ddg);
          }
        });
      });
    });

    // Collect all ΔG values
    Object.values(data.nodes).forEach(systemNodes => {
      systemNodes.forEach(node => {
        allDG.push(node.dg_exp);
        data.keys.forEach(key => {
          if (node.per_key[key]) {
            allDG.push(node.per_key[key].dg);
          }
        });
      });
    });

    // Calculate global domains with padding
    const ddgExtent = d3.extent(allDDG);
    const dgExtent = d3.extent(allDG);
    const ddgDomain = [ddgExtent[0] - 0.5, ddgExtent[1] + 0.5];
    const dgDomain = [dgExtent[0] - 0.5, dgExtent[1] + 0.5];

    // Two rows: DDG and DG
    scatterDiv.append('h4').text('ΔΔG vs Experiment');
    const ddgGrid = scatterDiv.append('div')
      .attr('class', 'scatter-row')
      .style('display', 'flex')
      .style('flex-wrap', 'wrap')
      .style('gap', `${layout.gap}px`);

    scatterDiv.append('h4').text('ΔG vs Experiment');
    const dgGrid = scatterDiv.append('div')
      .attr('class', 'scatter-row')
      .style('display', 'flex')
      .style('flex-wrap', 'wrap')
      .style('gap', `${layout.gap}px`);

    data.keys.forEach(key => {
      renderDDGScatter(ddgGrid.node(), systemId, key, ddgDomain, layout);
      renderDGScatter(dgGrid.node(), systemId, key, dgDomain, layout);
    });
  }

  function computeScatterLayout(container, keyCount) {
    const gap = 10;
    const minTile = 220;
    const maxTile = 440;
    const preferredTile = 320;
    const availableWidth = Math.max(minTile, getContainerWidth(container, minTile));

    const preferredCols = Math.max(1, Math.floor((availableWidth + gap) / (preferredTile + gap)));
    const columns = Math.max(1, Math.min(keyCount || 1, preferredCols));
    const rawTile = Math.floor((availableWidth - (columns - 1) * gap) / columns);
    const size = Math.max(minTile, Math.min(maxTile, rawTile));

    return {
      width: size,
      height: size,
      margin: { top: 30, right: 12, bottom: 36, left: 40 },
      pointRadius: Math.max(1.25, Math.min(2.5, Number((size * 0.008).toFixed(2)))),
      gap,
    };
  }

  function getContainerWidth(container, fallbackWidth) {
    const bounds = container.getBoundingClientRect();
    const width = Math.floor(bounds.width);
    return width > 0 ? width : fallbackWidth;
  }

  function computePairPlotLayout(container) {
    const gap = 20;
    const minWidth = 280;
    const maxWidth = 460;
    const preferredWidth = 360;
    const availableWidth = Math.max(minWidth, getContainerWidth(container, minWidth));
    const preferredColumns = Math.max(1, Math.floor((availableWidth + gap) / (preferredWidth + gap)));
    const columns = Math.max(1, Math.min(2, preferredColumns));
    const rawWidth = Math.floor((availableWidth - (columns - 1) * gap) / columns);
    const width = Math.max(minWidth, Math.min(maxWidth, rawWidth));

    return {
      width,
      height: width,
      gap,
      pointRadius: Math.max(1.5, Math.min(2.5, Number((width * 0.0065).toFixed(2)))),
    };
  }

  function renderDDGScatter(container, systemId, key, globalDomain, layout) {
    const width = layout.width;
    const height = layout.height;
    const margin = layout.margin;
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .style('width', `${width}px`)
      .style('height', `${height}px`)
      .attr('class', `scatter-ddg scatter-ddg-${key.replace(/\s+/g, '-')}`);

    const g = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // Title
    const labelHtml = renderLatex(key);
    const foreignObj = svg.append('foreignObject')
      .attr('x', margin.left)
      .attr('y', 5)
      .attr('width', plotWidth)
      .attr('height', 22);
    foreignObj.append('xhtml:div')
      .style('font-size', '10px')
      .style('font-weight', 'bold')
      .style('text-align', 'center')
      .html(labelHtml);

    // Data
    const edges = (data.edges[systemId] || []).filter(e => e.per_key[key]);
    if (edges.length === 0) {
      g.append('text')
        .attr('x', plotWidth / 2)
        .attr('y', plotHeight / 2)
        .attr('text-anchor', 'middle')
        .attr('font-size', '10px')
        .attr('fill', '#999')
        .text('No data');
      return;
    }

    // Use global domain for consistent axes across all plots
    const xScale = d3.scaleLinear().domain(globalDomain).range([0, plotWidth]);
    const yScale = d3.scaleLinear().domain(globalDomain).range([plotHeight, 0]);

    // Axes
    g.append('g')
      .attr('transform', `translate(0,${plotHeight})`)
      .call(d3.axisBottom(xScale).ticks(4).tickSize(3))
      .selectAll('text')
      .style('font-size', '8px');

    g.append('g')
      .call(d3.axisLeft(yScale).ticks(4).tickSize(3))
      .selectAll('text')
      .style('font-size', '8px');

    g.append('text')
      .attr('x', plotWidth / 2)
      .attr('y', plotHeight + 25)
      .attr('text-anchor', 'middle')
      .attr('font-size', '9px')
      .text('Exp (kcal/mol)');

    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -plotHeight / 2)
      .attr('y', -28)
      .attr('text-anchor', 'middle')
      .attr('font-size', '9px')
      .text('Calc (kcal/mol)');

    // Reference line
    g.append('line')
      .attr('x1', xScale(globalDomain[0]))
      .attr('x2', xScale(globalDomain[1]))
      .attr('y1', yScale(globalDomain[0]))
      .attr('y2', yScale(globalDomain[1]))
      .attr('stroke', '#ccc')
      .attr('stroke-width', 1)
      .attr('stroke-dasharray', '2,2');

    // Points
    const points = g.selectAll('.scatter-point')
      .data(edges)
      .join('circle')
      .attr('class', d => `scatter-point ddg-point ddg-point-${d.edge_id.replace(/\|/g, '-')}`)
      .attr('cx', d => xScale(d.ddg_exp))
      .attr('cy', d => yScale(d.per_key[key].ddg))
      .attr('r', layout.pointRadius)
      .attr('fill', colorScale(key))
      .attr('fill-opacity', 1)
      .attr('stroke', 'none')
      .attr('stroke-width', 0)
      .attr('cursor', 'pointer')
      .on('click', function(event, d) {
        handleEdgeSelection(systemId, d.edge_id, d.ligand_a, d.ligand_b);
      });

    // Background click to clear selection
    g.append('rect')
      .attr('width', plotWidth)
      .attr('height', plotHeight)
      .attr('fill', 'none')
      .attr('pointer-events', 'all')
      .lower()
      .on('click', function(event) {
        if (event.target === this) {
          clearSelection();
        }
      });

    updateSelectionOutlines();
  }

  function renderDGScatter(container, systemId, key, globalDomain, layout) {
    const width = layout.width;
    const height = layout.height;
    const margin = layout.margin;
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .style('width', `${width}px`)
      .style('height', `${height}px`)
      .attr('class', `scatter-dg scatter-dg-${key.replace(/\s+/g, '-')}`);

    const g = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // Title
    const labelHtml = renderLatex(key);
    const foreignObj = svg.append('foreignObject')
      .attr('x', margin.left)
      .attr('y', 5)
      .attr('width', plotWidth)
      .attr('height', 22);
    foreignObj.append('xhtml:div')
      .style('font-size', '10px')
      .style('font-weight', 'bold')
      .style('text-align', 'center')
      .html(labelHtml);

    // Data
    const nodes = (data.nodes[systemId] || []).filter(n => n.per_key[key]);
    if (nodes.length === 0) {
      g.append('text')
        .attr('x', plotWidth / 2)
        .attr('y', plotHeight / 2)
        .attr('text-anchor', 'middle')
        .attr('font-size', '10px')
        .attr('fill', '#999')
        .text('No data');
      return;
    }

    // Use global domain for consistent axes across all plots
    const xScale = d3.scaleLinear().domain(globalDomain).range([0, plotWidth]);
    const yScale = d3.scaleLinear().domain(globalDomain).range([plotHeight, 0]);

    // Axes
    g.append('g')
      .attr('transform', `translate(0,${plotHeight})`)
      .call(d3.axisBottom(xScale).ticks(4).tickSize(3))
      .selectAll('text')
      .style('font-size', '8px');

    g.append('g')
      .call(d3.axisLeft(yScale).ticks(4).tickSize(3))
      .selectAll('text')
      .style('font-size', '8px');

    g.append('text')
      .attr('x', plotWidth / 2)
      .attr('y', plotHeight + 25)
      .attr('text-anchor', 'middle')
      .attr('font-size', '9px')
      .text('Exp (kcal/mol)');

    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -plotHeight / 2)
      .attr('y', -28)
      .attr('text-anchor', 'middle')
      .attr('font-size', '9px')
      .text('Calc (kcal/mol)');

    // Reference line
    g.append('line')
      .attr('x1', xScale(globalDomain[0]))
      .attr('x2', xScale(globalDomain[1]))
      .attr('y1', yScale(globalDomain[0]))
      .attr('y2', yScale(globalDomain[1]))
      .attr('stroke', '#ccc')
      .attr('stroke-width', 1)
      .attr('stroke-dasharray', '2,2');

    // Points
    const points = g.selectAll('.scatter-point')
      .data(nodes)
      .join('circle')
      .attr('class', d => `scatter-point dg-point dg-point-${d.ligand.replace(/\s+/g, '-')}`)
      .attr('cx', d => xScale(d.dg_exp))
      .attr('cy', d => yScale(d.per_key[key].dg))
      .attr('r', layout.pointRadius)
      .attr('fill', colorScale(key))
      .attr('fill-opacity', 1)
      .attr('stroke', 'none')
      .attr('stroke-width', 0)
      .attr('cursor', 'pointer')
      .on('click', function(event, d) {
        handleLigandSelection(systemId, d.ligand);
      });

    g.append('rect')
      .attr('width', plotWidth)
      .attr('height', plotHeight)
      .attr('fill', 'none')
      .attr('pointer-events', 'all')
      .lower()
      .on('click', function(event) {
        if (event.target === this) {
          clearSelection();
        }
      });

    updateSelectionOutlines();
  }

  // ============================================================================
  // Selection Model
  // ============================================================================

  function handleEdgeSelection(systemId, edgeId, ligandA, ligandB) {
    if (selectedEdge && selectedEdge.edge_id === edgeId) {
      // Deselect
      clearSelection();
    } else {
      // Select
      selectedLigand = null;
      selectedEdge = { edge_id: edgeId, ligand_a: ligandA, ligand_b: ligandB, system_id: systemId };
      updateSelectionOutlines();
      updateStructureStrip(systemId, edgeId);
    }
  }

  function handleLigandSelection(systemId, ligand) {
    if (selectedLigand && selectedLigand.system_id === systemId && selectedLigand.ligand === ligand) {
      clearSelection();
      return;
    }

    selectedEdge = null;
    selectedLigand = { ligand, system_id: systemId };
    updateSelectionOutlines();
    updateLigandStructureStrip(systemId, ligand);
  }

  function clearSelection() {
    selectedEdge = null;
    selectedLigand = null;
    updateSelectionOutlines();
    clearStructureStrip();
  }

  function updateSelectionOutlines() {
    // Clear all outlines
    d3.selectAll('.scatter-point')
      .attr('stroke', 'none')
      .attr('stroke-width', 0);

    if (!selectedEdge && !selectedLigand) {
      updatePairPlotOutlines();
      return;
    }

    if (selectedLigand) {
      d3.selectAll('.ddg-point').each(function() {
        const point = d3.select(this);
        const pointData = point.datum();

        if (pointData.ligand_a === selectedLigand.ligand || pointData.ligand_b === selectedLigand.ligand) {
          point.attr('stroke', COLORS.purple)
               .attr('stroke-width', 4);
        }
      });

      d3.selectAll('.dg-point').each(function() {
        const point = d3.select(this);
        const pointData = point.datum();

        if (pointData.ligand === selectedLigand.ligand) {
          point.attr('stroke', COLORS.purple)
               .attr('stroke-width', 4);
        }
      });

      updatePairPlotOutlines();
      return;
    }

    // Apply outlines
    // DDG points
    d3.selectAll('.ddg-point').each(function() {
      const point = d3.select(this);
      const data = point.datum();

      if (data.edge_id === selectedEdge.edge_id) {
        // Purple outline for selected edge
        point.attr('stroke', COLORS.purple)
             .attr('stroke-width', 5);
      } else if (data.ligand_a === selectedEdge.ligand_a) {
        // Blue outline for ligand_a
        point.attr('stroke', COLORS.blue)
             .attr('stroke-width', 3.5);
      } else if (data.ligand_b === selectedEdge.ligand_b) {
        // Red outline for ligand_b
        point.attr('stroke', COLORS.red)
             .attr('stroke-width', 3.5);
      }
    });

    // DG points
    d3.selectAll('.dg-point').each(function() {
      const point = d3.select(this);
      const data = point.datum();

      if (data.ligand === selectedEdge.ligand_a) {
        point.attr('stroke', COLORS.blue)
             .attr('stroke-width', 3.5);
      } else if (data.ligand === selectedEdge.ligand_b) {
        point.attr('stroke', COLORS.red)
             .attr('stroke-width', 3.5);
      }
    });

    // Update pair plot if visible
    updatePairPlotOutlines();
  }

  // ============================================================================
  // Structure Strip
  // ============================================================================

  function updateStructureStrip(systemId, edgeId) {
    const strip = d3.select('#structure-strip');
    if (strip.empty()) return;

    strip.selectAll('*').remove();

    // Check if structures exist
    if (!data.structures || !data.structures[systemId] || !data.structures[systemId][edgeId]) {
      strip.append('p')
        .style('color', '#9ca3af')
        .style('font-style', 'italic')
        .text('No mapping available for this edge');
      return;
    }

    const structure = data.structures[systemId][edgeId];
    const edge = (data.edges[systemId] || []).find(e => e.edge_id === edgeId);
    if (!edge) return;

    // Edge info
    const info = strip.append('div')
      .style('margin-bottom', '10px');

    info.append('div')
      .style('font-weight', 'bold')
      .text(`Edge: ${edge.ligand_a} → ${edge.ligand_b}`);

    info.append('div')
      .style('font-size', '12px')
      .html(`ΔΔG<sub>exp</sub> = ${edge.ddg_exp.toFixed(2)} kcal/mol` +
            (edge.ddg_exp_unc ? ` ± ${edge.ddg_exp_unc.toFixed(2)}` : ''));

    // Structures
    const structContainer = strip.append('div')
      .style('display', 'flex')
      .style('gap', '10px')
      .style('justify-content', 'center')
      .style('align-items', 'center')
      .style('margin', '15px 0');

    const svgA = structContainer.append('div')
      .style('border', `3px solid ${COLORS.blue}`)
      .style('border-radius', '4px')
      .style('padding', '5px')
      .html(structure.svg_a);

    const arrow = structContainer.append('div')
      .style('font-size', '24px')
      .style('font-weight', 'bold')
      .text('→');

    const svgB = structContainer.append('div')
      .style('border', `3px solid ${COLORS.red}`)
      .style('border-radius', '4px')
      .style('padding', '5px')
      .html(structure.svg_b);

    // Caption (verbatim from plan)
    strip.append('p')
      .style('font-size', '11px')
      .style('color', '#6b7280')
      .style('font-style', 'italic')
      .style('margin-top', '10px')
      .text('Atom mapping from the benchmark ligand network (shared across all force fields — not this submission\'s mapping)');

    // Mapping stats
    strip.append('div')
      .style('font-size', '11px')
      .style('color', '#6b7280')
      .style('margin-top', '5px')
      .html(`Mapped: ${structure.n_mapped} atoms | ` +
            `Unique to ${edge.ligand_a}: ${structure.n_unique_a} | ` +
            `Unique to ${edge.ligand_b}: ${structure.n_unique_b}`);
  }

  function updateLigandStructureStrip(systemId, ligand) {
    const strip = d3.select('#structure-strip');
    if (strip.empty()) return;

    strip.selectAll('*').remove();

    const edges = data.edges[systemId] || [];
    const edge = edges.find(item => item.ligand_a === ligand || item.ligand_b === ligand);

    if (!edge || !data.structures || !data.structures[systemId] || !data.structures[systemId][edge.edge_id]) {
      strip.append('p')
        .style('color', '#9ca3af')
        .style('font-style', 'italic')
        .text('No structure available for this ligand');
      return;
    }

    const structure = data.structures[systemId][edge.edge_id];
    const isLigandA = edge.ligand_a === ligand;
    const ligandSvg = isLigandA ? structure.svg_a : structure.svg_b;

    strip.append('div')
      .style('font-weight', 'bold')
      .style('margin-bottom', '10px')
      .text(`Ligand: ${ligand}`);

    strip.append('div')
      .style('display', 'flex')
      .style('justify-content', 'center')
      .style('margin', '15px 0')
      .append('div')
      .style('border', `3px solid ${COLORS.purple}`)
      .style('border-radius', '4px')
      .style('padding', '5px')
      .html(ligandSvg);

    strip.append('p')
      .style('font-size', '11px')
      .style('color', '#6b7280')
      .style('font-style', 'italic')
      .style('margin-top', '10px')
      .text('Single-ligand depiction taken from one available mapped edge for the selected ligand.');
  }

  function clearStructureStrip() {
    const strip = d3.select('#structure-strip');
    if (strip.empty()) return;

    strip.selectAll('*').remove();
    strip.append('div')
      .attr('class', 'structure-placeholder')
      .style('color', '#9ca3af')
      .style('font-style', 'italic')
      .text('Click a scatter point to view 2D structures');
  }

  // ============================================================================
  // Pair Plot Control
  // ============================================================================

  function renderPairPlotControl(container, systemId) {
    const div = d3.select(container);

    // Controls
    const controls = div.append('div')
      .style('margin-bottom', '15px')
      .style('display', 'flex')
      .style('gap', '15px')
      .style('align-items', 'center');

    controls.append('label').text('X-axis:');
    const key1Select = controls.append('select')
      .attr('id', 'pair-key1')
      .style('padding', '5px');

    controls.append('label').text('Y-axis:');
    const key2Select = controls.append('select')
      .attr('id', 'pair-key2')
      .style('padding', '5px');

    controls.append('button')
      .text('Plot')
      .style('padding', '5px 15px')
      .style('cursor', 'pointer')
      .on('click', () => {
        const key1 = key1Select.node().value;
        const key2 = key2Select.node().value;
        if (key1 && key2 && key1 !== key2) {
          renderPairPlot(systemId, key1, key2);
        }
      });

    // Populate selects
    data.keys.forEach((key, i) => {
      key1Select.append('option').attr('value', key).text(key);
      key2Select.append('option').attr('value', key).text(key);
    });

    if (data.keys.length >= 2) {
      key1Select.node().selectedIndex = 0;
      key2Select.node().selectedIndex = 1;
    }

    // Pair plot container
    div.append('div').attr('id', 'pair-plot-charts');
  }

  function renderPairPlot(systemId, key1, key2) {
    const container = d3.select('#pair-plot-charts');
    container.selectAll('*').remove();
    const layout = computePairPlotLayout(container.node());

    const plotDiv = container.append('div')
      .style('display', 'flex')
      .style('flex-wrap', 'wrap')
      .style('gap', `${layout.gap}px`);

    renderDDGPairPlot(plotDiv.node(), systemId, key1, key2, layout);
    renderDGPairPlot(plotDiv.node(), systemId, key1, key2, layout);

    updatePairPlotOutlines();
  }

  function renderDDGPairPlot(container, systemId, key1, key2, layout) {
    const width = layout.width;
    const height = layout.height;
    const margin = { top: 30, right: 10, bottom: 50, left: 50 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height);

    const g = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    svg.append('text')
      .attr('x', width / 2)
      .attr('y', 15)
      .attr('text-anchor', 'middle')
      .attr('font-weight', 'bold')
      .text('ΔΔG Pair Plot');

    // Data: shared edges
    const edges = (data.edges[systemId] || [])
      .filter(e => e.per_key[key1] && e.per_key[key2]);

    const excluded = (data.edges[systemId] || []).length - edges.length;

    if (edges.length === 0) {
      g.append('text')
        .attr('x', plotWidth / 2)
        .attr('y', plotHeight / 2)
        .attr('text-anchor', 'middle')
        .attr('font-size', '10px')
        .attr('fill', '#999')
        .text('No shared edges');
      return;
    }

    const allVals = edges.map(e => [e.per_key[key1].ddg, e.per_key[key2].ddg]).flat();
    const extent = d3.extent(allVals);
    const domain = [extent[0] - 0.5, extent[1] + 0.5];

    const xScale = d3.scaleLinear().domain(domain).range([0, plotWidth]);
    const yScale = d3.scaleLinear().domain(domain).range([plotHeight, 0]);

    // Axes
    g.append('g')
      .attr('transform', `translate(0,${plotHeight})`)
      .call(d3.axisBottom(xScale).ticks(5));

    g.append('g')
      .call(d3.axisLeft(yScale).ticks(5));

    g.append('text')
      .attr('x', plotWidth / 2)
      .attr('y', plotHeight + 40)
      .attr('text-anchor', 'middle')
      .attr('font-size', '11px')
      .text(`${key1} (kcal/mol)`);

    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -plotHeight / 2)
      .attr('y', -35)
      .attr('text-anchor', 'middle')
      .attr('font-size', '11px')
      .text(`${key2} (kcal/mol)`);

    // Reference line
    g.append('line')
      .attr('x1', xScale(domain[0]))
      .attr('x2', xScale(domain[1]))
      .attr('y1', yScale(domain[0]))
      .attr('y2', yScale(domain[1]))
      .attr('stroke', '#ccc')
      .attr('stroke-width', 1)
      .attr('stroke-dasharray', '2,2');

    // Points
    g.selectAll('.pair-point')
      .data(edges)
      .join('circle')
      .attr('class', d => `scatter-point pair-ddg-point pair-ddg-point-${d.edge_id.replace(/\|/g, '-')}`)
      .attr('cx', d => xScale(d.per_key[key1].ddg))
      .attr('cy', d => yScale(d.per_key[key2].ddg))
      .attr('r', layout.pointRadius)
      .attr('fill', COLORS.gray)
      .attr('fill-opacity', 1)
      .attr('stroke', 'none')
      .attr('stroke-width', 0)
      .attr('cursor', 'pointer')
      .on('click', function(event, d) {
        handleEdgeSelection(systemId, d.edge_id, d.ligand_a, d.ligand_b);
      });

    if (excluded > 0) {
      svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 5)
        .attr('text-anchor', 'middle')
        .attr('font-size', '10px')
        .attr('fill', '#ef4444')
        .text(`${excluded} edge${excluded !== 1 ? 's' : ''} excluded (not in both keys)`);
    }
  }

  function renderDGPairPlot(container, systemId, key1, key2, layout) {
    const width = layout.width;
    const height = layout.height;
    const margin = { top: 30, right: 10, bottom: 50, left: 50 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height);

    const g = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    svg.append('text')
      .attr('x', width / 2)
      .attr('y', 15)
      .attr('text-anchor', 'middle')
      .attr('font-weight', 'bold')
      .text('ΔG Pair Plot');

    // Data: shared nodes
    const nodes = (data.nodes[systemId] || [])
      .filter(n => n.per_key[key1] && n.per_key[key2]);

    const excluded = (data.nodes[systemId] || []).length - nodes.length;

    if (nodes.length === 0) {
      g.append('text')
        .attr('x', plotWidth / 2)
        .attr('y', plotHeight / 2)
        .attr('text-anchor', 'middle')
        .attr('font-size', '10px')
        .attr('fill', '#999')
        .text('No shared nodes');
      return;
    }

    const allVals = nodes.map(n => [n.per_key[key1].dg, n.per_key[key2].dg]).flat();
    const extent = d3.extent(allVals);
    const domain = [extent[0] - 0.5, extent[1] + 0.5];

    const xScale = d3.scaleLinear().domain(domain).range([0, plotWidth]);
    const yScale = d3.scaleLinear().domain(domain).range([plotHeight, 0]);

    // Axes
    g.append('g')
      .attr('transform', `translate(0,${plotHeight})`)
      .call(d3.axisBottom(xScale).ticks(5));

    g.append('g')
      .call(d3.axisLeft(yScale).ticks(5));

    g.append('text')
      .attr('x', plotWidth / 2)
      .attr('y', plotHeight + 40)
      .attr('text-anchor', 'middle')
      .attr('font-size', '11px')
      .text(`${key1} (kcal/mol)`);

    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', -plotHeight / 2)
      .attr('y', -35)
      .attr('text-anchor', 'middle')
      .attr('font-size', '11px')
      .text(`${key2} (kcal/mol)`);

    // Reference line
    g.append('line')
      .attr('x1', xScale(domain[0]))
      .attr('x2', xScale(domain[1]))
      .attr('y1', yScale(domain[0]))
      .attr('y2', yScale(domain[1]))
      .attr('stroke', '#ccc')
      .attr('stroke-width', 1)
      .attr('stroke-dasharray', '2,2');

    // Points
    g.selectAll('.pair-point')
      .data(nodes)
      .join('circle')
      .attr('class', d => `scatter-point pair-dg-point pair-dg-point-${d.ligand.replace(/\s+/g, '-')}`)
      .attr('cx', d => xScale(d.per_key[key1].dg))
      .attr('cy', d => yScale(d.per_key[key2].dg))
      .attr('r', layout.pointRadius)
      .attr('fill', COLORS.gray)
      .attr('fill-opacity', 1)
      .attr('stroke', 'none')
      .attr('stroke-width', 0)
      .attr('cursor', 'pointer')
      .on('click', function(event, d) {
        handleLigandSelection(systemId, d.ligand);
      });

    if (excluded > 0) {
      svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 5)
        .attr('text-anchor', 'middle')
        .attr('font-size', '10px')
        .attr('fill', '#ef4444')
        .text(`${excluded} node${excluded !== 1 ? 's' : ''} excluded (not in both keys)`);
    }
  }

  function updatePairPlotOutlines() {
    d3.selectAll('.pair-ddg-point, .pair-dg-point')
      .attr('stroke', 'none')
      .attr('stroke-width', 0);

    if (!selectedEdge && !selectedLigand) return;

    if (selectedLigand) {
      d3.selectAll('.pair-ddg-point').each(function() {
        const point = d3.select(this);
        const pointData = point.datum();

        if (pointData.ligand_a === selectedLigand.ligand || pointData.ligand_b === selectedLigand.ligand) {
          point.attr('stroke', COLORS.purple).attr('stroke-width', 4);
        }
      });

      d3.selectAll('.pair-dg-point').each(function() {
        const point = d3.select(this);
        const pointData = point.datum();

        if (pointData.ligand === selectedLigand.ligand) {
          point.attr('stroke', COLORS.purple).attr('stroke-width', 4);
        }
      });

      return;
    }

    // Update pair plot DDG points
    d3.selectAll('.pair-ddg-point').each(function() {
      const point = d3.select(this);
      const pointData = point.datum();

      if (pointData.edge_id === selectedEdge.edge_id) {
        point.attr('stroke', COLORS.purple).attr('stroke-width', 5);
      } else if (pointData.ligand_a === selectedEdge.ligand_a) {
        point.attr('stroke', COLORS.blue).attr('stroke-width', 3.5);
      } else if (pointData.ligand_b === selectedEdge.ligand_b) {
        point.attr('stroke', COLORS.red).attr('stroke-width', 3.5);
      }
    });

    // Update pair plot DG points
    d3.selectAll('.pair-dg-point').each(function() {
      const point = d3.select(this);
      const pointData = point.datum();

      if (pointData.ligand === selectedEdge.ligand_a) {
        point.attr('stroke', COLORS.blue).attr('stroke-width', 3.5);
      } else if (pointData.ligand === selectedEdge.ligand_b) {
        point.attr('stroke', COLORS.red).attr('stroke-width', 3.5);
      }
    });
  }

  // ============================================================================
  // Initialization
  // ============================================================================

  function init() {
    try {
      // Load data
      data = loadData();

      // Setup color scale
      colorScale = setupColorScale(data.keys);

      // Render main panel
      const root = document.getElementById('widget-root');
      renderMainPanel(root);

      console.log('Widget initialized successfully');
    } catch (error) {
      console.error('Widget initialization failed:', error);
      const root = document.getElementById('widget-root');
      root.innerHTML = `<div style="padding: 20px; color: red;">
        <h2>Error loading widget</h2>
        <p>${error.message}</p>
      </div>`;
    }
  }

  // Run on DOM ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

})();
