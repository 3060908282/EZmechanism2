#!/bin/bash
cd /home/z/my-project/mini-services/prediction-service
while true; do
    echo "[$(date)] Starting prediction service..."
    python3 -u index.py 2>&1
    echo "[$(date)] Service exited with code $?, restarting in 3s..."
    sleep 3
done
