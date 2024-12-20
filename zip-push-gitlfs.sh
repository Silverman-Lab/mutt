#!/bin/bash

# Repository Remote URL
REMOTE_URL="git@github.com:Silverman-Lab/data_repository.git"

# Ensure Git and Git LFS are installed
if ! command -v git &> /dev/null || ! command -v git-lfs &> /dev/null; then
    echo "Error: Git and/or Git LFS is not installed. Please install them and try again."
    exit 1
fi

# Ensure the current directory is a Git repository
if [ ! -d ".git" ]; then
    echo "Error: Current directory is not a Git repository. Please cd into the repository's root and try again."
    exit 1
fi

# Verify or Set the Remote Repository
CURRENT_REMOTE=$(git remote get-url origin 2>/dev/null || echo "")
if [ "$CURRENT_REMOTE" != "$REMOTE_URL" ]; then
    echo "Setting the correct remote repository..."
    git remote remove origin 2>/dev/null || echo "No existing remote to remove."
    git remote add origin "$REMOTE_URL"
fi

# Confirm the remote URL
echo "Remote repository is set to:"
git remote -v

# Update the repository while preserving local changes
echo "Updating the repository..."
git stash save "Backup before pull" &> /dev/null || echo "No changes to stash."
git pull origin main || echo "Error pulling from remote repository. Please check the connection."
git stash pop &> /dev/null || echo "No stash to apply."

# Process repository files
echo "Processing repository files..."
BASE_DIR=$(pwd)

# Find all subdirectories excluding .git
find . -type d -not -path "./.git*" | while read -r folder; do
    # Skip directories containing .gitattributes
    if [ -f "$folder/.gitattributes" ]; then
        echo "Skipping folder $folder as .gitattributes is present."
        continue
    fi

    echo "Processing folder: $folder"
    cd "$folder" || continue

    # Zip files excluding README variants, parse.R, *.zip, and *.gz
    for file in *; do
        if [ -f "$file" ]; then
            case "$file" in
                # Matches README, readme, Readme, etc. and any extensions following it (e.g. README.md, readme.txt)
                [Rr][Ee][Aa][Dd][Mm][Ee]*|parse.R|*.zip|*.gz)
                    # Skip these files
                    ;;
                *)
                    echo "Zipping $file in $folder"
                    zip "${file}.zip" "$file"
                    rm -f "$file"
                    ;;
            esac
        fi
    done

    cd "$BASE_DIR" || exit
done

# Track .zip files with Git LFS
echo "Tracking .zip files with Git LFS..."
git lfs track "*.zip"

# Stage changes
echo "Staging changes..."
git add .

# Commit changes if there are any
if git diff-index --quiet HEAD; then
    echo "No changes to commit."
else
    git commit -m "Zip files and track with Git LFS"
fi

# Push changes to the remote repository
echo "Pushing changes to the repository..."
git push origin main || echo "Error pushing to remote repository. Please check your connection."

