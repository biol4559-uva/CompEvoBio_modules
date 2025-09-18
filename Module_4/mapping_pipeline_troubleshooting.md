# **mapping pipeline troubleshooting**

## Overview
For some folks were having difficulties getting their samples to map. In some cases, this has to do with path names. In other cases, this has to do with special hidden characters that Windows uses. You'll need to delete those characters form your "small_file.csv".

You can see if your output got flagged as not working if it ended up in this folder `/standard/BerglandTeach/badOutput`.

If so, try running the following line on your "small_file.csv". You'll need to navigate into the directory where this file exists first.

```cat small_file.csv | sed $'s/[^[:print:]\t]//g' | tr -d '^M' > small_file_clean.csv```

One you run this, try
```cat -A small_file.csv```
```cat -A small_file_clean.csv```

That shoudl help you determine if your new file is cleaned up or not. If not, flag us down to chat.
