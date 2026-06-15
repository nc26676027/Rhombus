#!/bin/bash
# Run all Rhombus tests and exit non-zero on any failure.
# Usage: bash scripts/run-tests.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN_DIR="$PROJECT_DIR/build/bin"

PASS=0
FAIL=0

run_test() {
    local name="$1"
    local bin="$BIN_DIR/$name"
    echo "--- Running: $name ---"
    if [ ! -f "$bin" ]; then
        echo "SKIP: $bin not found"
        return
    fi
    if "$bin" 2>&1 | tee /dev/stderr | grep -q "max_error\|max error\|Passed"; then
        # Check that max_error values are within acceptable range (<=2 for PP-MVM)
        local max_err
        max_err=$("$bin" 2>&1 | grep -oP 'max_error\s*=\s*\K\d+|max error:\s*\K\d+' | sort -n | tail -1)
        if [ -n "$max_err" ] && [ "$max_err" -le 2 ]; then
            echo "PASS: $name (max_error=$max_err)"
            PASS=$((PASS + 1))
        else
            echo "FAIL: $name (max_error=$max_err)"
            FAIL=$((FAIL + 1))
        fi
    else
        echo "FAIL: $name (no output or crash)"
        FAIL=$((FAIL + 1))
    fi
    echo ""
}

# Core tests
run_test "matvec"
run_test "matmul"

# Face module tests (Phase 1+)
[ -f "$BIN_DIR/test_hydia" ] && run_test "test_hydia"
[ -f "$BIN_DIR/test_threshold" ] && run_test "test_threshold"
[ -f "$BIN_DIR/test_tagselect" ] && run_test "test_tagselect"

echo "=== Summary ==="
echo "PASS: $PASS  FAIL: $FAIL"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
