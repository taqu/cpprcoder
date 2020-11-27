# cpprcoder
An implementaion of adaptive range coder for C++

# Benchmarks
## System
 Windows 10 Pro Ver.2004(Build 19041.508) 64bit  
 Visual Studio Professional 2019 Version 16.7.5  
 CPU: Intel Core i7-8700  
 Memory: 32 GB (DDR4)

## Test data
- [The Canterbury Corpus](http://corpus.canterbury.ac.nz/index.html)

## Result table

### Range Coder

|Name|Ratio|Compression (micro seconds)|Decompression (micro seconds)|
|:---|:---|:---|:---|
|alice29.txt|0.574532|1674|5964|
|asyoulik.txt|0.605293|1416|4944|
|cp.html|0.674836|276|957|
|fields.c|0.672646|129|437|
|grammar.lsp|0.718893|43|147|
|kennedy.xls|0.452938|10656|39724|
|lcet10.txt|0.585129|4698|16673|
|plrabn12.txt|0.567788|5370|18808|
|ptt5|0.157010|4946|19490|
|sum|0.679759|481|1533|
|xargs.1|0.735510|48|167|

### Adaptive Range Coder

|Name|Ratio|Compression (micro seconds)|Decompression (micro seconds)|
|:---|:---|:---|:---|
|alice29.txt|0.573000|3224|4966|
|asyoulik.txt|0.603400|2716|4360|
|cp.html|0.662480|569|857|
|fields.c|0.642511|244|367|
|grammar.lsp|0.619457|83|118|
|kennedy.xls|0.447426|18369|22158|
|lcet10.txt|0.584625|9151|14856|
|plrabn12.txt|0.567367|11152|16143|
|ptt5|0.152158|7112|9837|
|sum|0.670450|814|1164|
|xargs.1|0.648924|98|149|

### SLZ4
|Name|Ratio|Compression (micro seconds)|Decompression (micro seconds)|
|:---|:---|:---|:---|
|alice29.txt|0.590510|1202|406|
|asyoulik.txt|0.623499|1016|338|
|cp.html|0.494452|184|45|
|fields.c|0.479283|73|22|
|grammar.lsp|0.527546|30|6|
|kennedy.xls|0.362974|4801|1208|
|lcet10.txt|0.552485|3124|1078|
|plrabn12.txt|0.654747|4174|1360|
|sum|0.501464|288|67|
|xargs.1|0.637568|40|8|

### ZLib
|Name|Ratio|Compression (micro seconds)|Decompression (micro seconds)|
|:---|:---|:---|:---|
|alice29.txt|0.357712|12499|1598|
|asyoulik.txt|0.390617|11068|1325|
|cp.html|0.323578|1094|262|
|fields.c|0.280000|734|124|
|grammar.lsp|0.328406|241|54|
|kennedy.xls|0.198100|62349|8008|
|lcet10.txt|0.339549|39025|4607|
|plrabn12.txt|0.405223|55098|5270|
|ptt5|0.110022|20496|3429|
|sum|0.339697|2839|389|
|xargs.1|0.410693|245|65|

### LZ4
|Name|Ratio|Compression (micro seconds)|Decompression (micro seconds)|
|:---|:---|:---|:---|
|alice29.txt|0.583205|509|189|
|asyoulik.txt|0.636313|447|102|
|cp.html|0.483884|50|16|
|fields.c|0.467713|21|7|
|grammar.lsp|0.513840|7|1|
|kennedy.xls|0.363881|1371|532|
|lcet10.txt|0.546481|1034|345|
|plrabn12.txt|0.675691|1280|271|
|ptt5|0.169295|411|306|
|sum|0.491946|113|22|
|xargs.1|0.628815|10|1|

# Building
First, make sure you clone this repository.  
Then, deflate cantrbry.tar.bz2. 
If you work on Windows, you're going to need to deflate zlib.tar.bz2 and lz4.tar.bz2.

```
git clone https://github.com/taqu/cpprcoder.git
cd cpprcoder/test
tar jxvf cantrbry.tar.bz2
```

Next,
```
mkdir build
pushd build
cmake -DUSE_SIMD=1 ..
```

# License
This software is distributed under two licenses 'The MIT License' or 'Public Domain', choose whichever you like.
