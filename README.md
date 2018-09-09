# cpparcoder
An implementaion of adaptive range coder for C++

# Benchmarks
## System
 Windows 10 Pro Ver.1803(Build 17134.228) 64bit  
 Visual Studio 2017 Ver 15.5.2  
 CPU: Intel Core i7-6700T  
 Memory: 16GB (DDR4)

## Test data
- [The Canterbury Corpus](http://corpus.canterbury.ac.nz/index.html)
## Result table
### Adaptive Range Coder
|Name|Ratio|Compression (nano seconds)|Decompression (nano seconds)|
|:---|:---|:---|:---|
|alice29.txt|0.573000|5587|7150|
|asyoulik.txt|0.603400|3769|5472|
|cp.html|0.662480|785|1314|
|fields.c|0.642511|476|582|
|grammar.lsp|0.619457|110|168|
|kennedy.xls|0.447426|25035|29346|
|lcet10.txt|0.584625|13651|18465|
|plrabn12.txt|0.567367|18652|21234|
|ptt5|0.152158|10315|13233|
|sum|0.670450|1104|1530|
|xargs.1|0.648924|129|188|

### zlib
|Name|Ratio|Compression (nano seconds)|Decompression (nano seconds)|
|:---|:---|:---|:---|
|alice29.txt|0.357712|7935|8979|
|asyoulik.txt|0.390617|6997|7856|
|cp.html|0.323578|708|871|
|fields.c|0.280000|324|401|
|grammar.lsp|0.328406|147|183|
|kennedy.xls|0.198100|33404|39149|
|lcet10.txt|0.339549|21247|23925|
|plrabn12.txt|0.405223|32987|36529|
|ptt5|0.110022|9962|11931|
|sum|0.339697|1669|1940|
|xargs.1|0.410693|167|209|

### LZ4
|Name|Ratio|Compression (nano seconds)|Decompression (nano seconds)|
|:---|:---|:---|:---|
|alice29.txt|0.583205|477|564|
|asyoulik.txt|0.636313|365|433|
|cp.html|0.483884|61|76|
|fields.c|0.467713|33|44|
|grammar.lsp|0.513840|12|14|
|kennedy.xls|0.363881|1934|2454|
|lcet10.txt|0.546481|1235|1466|
|plrabn12.txt|0.675691|1457|1715|
|ptt5|0.169295|513|780|
|sum|0.491946|88|112|
|xargs.1|0.628815|13|15|

# License
This is free and unencumbered software released into the public domain.
