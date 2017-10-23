import sys
def progress(kSample,NSAMPLE,preamble=''):
        a = kSample*100/NSAMPLE
        b = (NSAMPLE-1-kSample)*100/NSAMPLE
        sys.stdout.write('\r'+preamble+'|'+'~'*a+' '*b+'|')
