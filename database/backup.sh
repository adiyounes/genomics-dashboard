#!/bin/bash
# Run this weekly: bash database/backup.sh
pg_dump -h localhost -U genomics_user genomics_db > database/backups/genomics_db_$(date +%Y%m%d).sql
echo "✅ Backup saved"
