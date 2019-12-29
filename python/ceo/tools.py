import sys
def progress(kSample,NSAMPLE,preamble=''):
    a = kSample*100/NSAMPLE
    b = (NSAMPLE-1-kSample)*100/NSAMPLE
    sys.stdout.write('\r'+preamble+'|'+'~'*a+' '*b+'|')

def ascupy(ceo_cu_array):
    from cupy.cuda import UnownedMemory, MemoryPointer
    from cupy import ndarray
    ptr = UnownedMemory(ceo_cu_array.dev_ptr,
                        ceo_cu_array.nbytes,
                        ceo_cu_array)
    return ndarray(shape=ceo_cu_array.shape,
                      dtype=ceo_cu_array.type,
                      memptr=MemoryPointer(ptr,0))
