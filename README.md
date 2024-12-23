#Analysis

For the first step i used MACS2 in my terminal to find Peaks, this is the code i used :

For file in `Task/BED/*.tagAlign` do:
```bash
macs2 callpeak -t "$file" -f AUTO -g hs -n "$(basename "$file")" --outdir /mnt/c/Users/Administrator/Desktop/Task/Peaks --nomodel --extsize 137
```

Then i created the unionpeak using bedtools:
```bash
cat /mnt/c/Users/Administrator/Desktop/Task/Peaks/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > /mnt/c/Users/Administrator/Desktop/Task/Peaks/union_peaks.bed
```

The rest was done in R. You will find two scripts, one named Training and one named Final.
The Trainig script is the script i wrote for the first time without using any loop to get a clear picture of what i was doing.
I decided to upload it to show you the logic i followed and my progress
The Final script is the updated one
