#!/bin/bash
# Regenerate and open the comparison widget

echo "Regenerating comparison widget..."

REFERENCE_MODEL="Sage 2.3.0 (ff14SB) TIP3P"

# Run the export script
micromamba run -n openfe python ../../src/export_comparison_widget_data.py \
  --config comparison_widget_config.yaml \
  --output-dir . \
  --num-bootstraps 5000 \
  --ecdf-bootstraps 0 \
  --seed 42 \
  --reference-model "$REFERENCE_MODEL"

#  --no-structures \

echo ""
echo "Widget generated successfully!"

# Open the widget in the default browser
WIDGET_PATH="comparison_widget.html"

if [ -f "$WIDGET_PATH" ]; then
  echo "Opening $WIDGET_PATH in browser..."

  # Detect OS and open accordingly
  if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    open "$WIDGET_PATH"
  elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux
    xdg-open "$WIDGET_PATH" 2>/dev/null || echo "Please open $WIDGET_PATH manually"
  else
    echo "Please open $WIDGET_PATH manually"
  fi
else
  echo "Error: Widget file not found at $WIDGET_PATH"
  exit 1
fi
