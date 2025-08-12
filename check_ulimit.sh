#!/bin/bash

# Script to check and optionally increase file descriptor limits
# This addresses the "Too many open files" error in BindCraft

echo "BindCraft File Descriptor Limit Checker"
echo "========================================"

# Check current limits
echo "Current file descriptor limits:"
echo "Soft limit: $(ulimit -n)"
echo "Hard limit: $(ulimit -Hn)"

# Check system-wide limits
if [ -f /proc/sys/fs/file-max ]; then
    echo "System-wide limit: $(cat /proc/sys/fs/file-max)"
fi

# Recommended minimum for BindCraft
RECOMMENDED_LIMIT=65536

current_soft=$(ulimit -n)
if [ "$current_soft" -lt "$RECOMMENDED_LIMIT" ]; then
    echo ""
    echo "WARNING: Current soft limit ($current_soft) is below recommended ($RECOMMENDED_LIMIT)"
    echo ""
    echo "To fix the 'Too many open files' error, you can:"
    echo ""
    echo "1. TEMPORARY FIX (current session only):"
    echo "   ulimit -n $RECOMMENDED_LIMIT"
    echo ""
    echo "2. PERMANENT FIX (add to ~/.bashrc or ~/.zshrc):"
    echo "   echo 'ulimit -n $RECOMMENDED_LIMIT' >> ~/.bashrc"
    echo "   source ~/.bashrc"
    echo ""
    echo "3. SYSTEM-WIDE FIX (requires root, add to /etc/security/limits.conf):"
    echo "   * soft nofile $RECOMMENDED_LIMIT"
    echo "   * hard nofile $RECOMMENDED_LIMIT"
    echo ""
    
    read -p "Would you like to set the temporary limit now? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        ulimit -n $RECOMMENDED_LIMIT
        echo "File descriptor limit temporarily set to $(ulimit -n)"
        echo "You can now run BindCraft in this terminal session."
    fi
else
    echo ""
    echo "✓ File descriptor limit looks good for BindCraft!"
fi

echo ""
echo "For more information about this issue, see:"
echo "https://github.com/your-repo/BindCraft/issues/file-descriptor-leak"
