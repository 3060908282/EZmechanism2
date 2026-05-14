#!/bin/bash
# Start all services: prediction (3003), ketcher (3004), next.js (3000)
set +e

# Start prediction service in background
(
  cd /home/z/my-project/mini-services/prediction-service
  while true; do
    python3 -u index.py 2>&1
    echo "[prediction] Restarting in 3s..."
    sleep 3
  done
) &
echo "[dev-all] Prediction service started (PID: $!)"

# Start ketcher service in background
(
  cd /home/z/my-project/mini-services/ketcher-service
  while true; do
    node server.js 2>&1
    echo "[ketcher] Restarting in 3s..."
    sleep 3
  done
) &
echo "[dev-all] Ketcher service started (PID: $!)"

# Start Next.js in foreground (this keeps the script alive)
cd /home/z/my-project
rm -rf .next
exec node node_modules/.bin/next dev -p 3000
