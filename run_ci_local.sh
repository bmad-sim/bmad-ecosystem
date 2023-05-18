#!/bin/bash
act workflow_dispatch --container-architecture linux/amd64 -W .github/workflows/ci.yml
