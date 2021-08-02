
import io
import codecs

def detect_by_bom(path,default):
    with open(path, 'rb') as f:
        raw = f.read(32)    #will read less if the file is smaller
        for enc,boms in \
                ('utf-8-sig',(codecs.BOM_UTF8,)),\
                ('utf-16',(codecs.BOM_UTF16_LE,codecs.BOM_UTF16_BE)),\
                ('utf-32',(codecs.BOM_UTF32_LE,codecs.BOM_UTF32_BE)):
            if any(raw.startswith(bom) for bom in boms): return enc

    return default

def cat(file):
    enc=detect_by_bom(file,'utf8')

    with io.open(file, encoding=enc, errors='ignore') as f:
        return f.read()


