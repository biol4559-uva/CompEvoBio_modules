# **mapping pipeline troubleshooting**

## Overview
For some folks were having difficulties getting their samples to map. In some cases, this has to do with path names. In other cases, this has to do with special hidden characters that Windows uses. You'll need to delete those characters form your "small_file.csv".

You can see if your output got flagged as not working if it ended up in this folder `/standard/BerglandTeach/badOutput`.

If so, try running the following line, but using the path to your small file:

```cat small_file.csv | $'s/[^[:print:]\t]//g' | tr -d '^M'```
