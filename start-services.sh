#!/bin/bash
# Persistent service launcher for all EZmechanism services
# Each service runs in its own auto-restart loop
set +e

PREDICTION_DIR="/home/z/my-project/mini-services/prediction-service"
KETCHER_DIR="/home/z/my-project/mini-services/ketcher-service"
PROJECT_DIR="/home/z/my-project"

echo "=== Starting EZmechanism Services ==="

# Start prediction service (Flask on port 3003)
(
  cd "$PREDICTION_DIR"
  while true; do
    echo "[$(date '+%H:%M:%S')] Starting prediction service on port 3003..."
    /home/z/.venv/bin/python3 -u index.py 2>&1
    echo "[$(date '+%H:%M:%S')] Prediction service exited, restarting in 3s..."
    sleep 3
  done
) &
PREDICTION_PID=$!
echo "Prediction service PID: $PREDICTION_PID"

# Start ketcher service (Node on port 3004)
(
  cd "$KETCHER_DIR"
  while true; do
    echo "[$(date '+%H:%M:%S')] Starting ketcher service on port 3004..."
    node server.js 2>&1
    echo "[$(date '+%H:%M:%S')] Ketcher service exited, restarting in 3s..."
    sleep 3
  done
) &
KETCHER_PID=$!
echo "Ketcher service PID: $KETCHER_PID"

# Start Next.js dev server (port 3000)
# Do NOT delete .next — that causes 8s+ recompilation on every restart
(
  cd "$PROJECT_DIR"
  while true; do
    echo "[$(date '+%H:%M:%S')] Starting Next.js dev server on port 3000..."
    node node_modules/.bin/next dev -p 3000 2>&1
    echo "[$(date '+%H:%M:%S')] Next.js exited, restarting in 5s..."
    sleep 5
  done
) &
NEXTJS_PID=$!
echo "Next.js PID: $NEXTJS_PID"

echo ""
echo "=== All services launched ==="
echo "  Prediction (3003): PID $PREDICTION_PID"
echo "  Ketcher    (3004): PID $KETCHER_PID"
echo "  Next.js    (3000): PID $NEXTJS_PID"
echo ""
echo "Services will auto-restart if they crash."

# Wait for all background processes
wait
