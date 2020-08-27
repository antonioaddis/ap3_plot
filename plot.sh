#1 tstart
#2 tstop
#3 first vertical line
#4 second vertical line

echo "~/projects/IGRJ17354-3255/2020/Table_3/dir_100_60.0_$1.0_$2.0"

#python plot.py --agile ~/projects/IGRJ17354-3255/2020/LS/AGILE/LSclean4/merged/RES1_TBS3600_R2_DQ1_EB100-50000.ap.ap4 --fermi ~/projects/IGRJ17354-3255/2020/LS/FERMI/FER1_TBS3600_R2_BARC0_EB100-10000.ap.ap4 --tstart $1 --tstop $2 --path_offaxis "/Users/bulgarelli/projects/IGRJ17354-3255/2020/Table_3/dir_100_60.0_$1.0_$2.0" --line1 $3 --line2 $4



python plot.py --agile ~/projects/IGRJ17354-3255/2020/LS/AGILE/LSclean4/merged/RES1_TBS3600_R2_DQ1_EB100-50000.ap.ap4 --fermi /Users/bulgarelli/docker/fermibootle/FER1_TBS300_R2_BARC0_EB100-10000_MJD55598_55601.ap.ap4 --tstart $1 --tstop $2 --path_offaxis "/Users/bulgarelli/projects/IGRJ17354-3255/2020/Table_3/dir_100_60.0_$1.0_$2.0" --line1 $3 --line2 $4


