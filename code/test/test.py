import ctypes  
ll = ctypes.cdll.LoadLibrary   
lib = ll("./libfun.so")

#lib = ll("./libpycall.so")    
lib_mine = ll("./libTestMine.so")
#lib.foo(1, 3) 
s = 'in.root'
p = ctypes.c_wchar_p(s)
b = lib_mine.string_conventor(s.encode('ascii'),len(s))
lib.reasdd(b.encode('ascii'),len(b))
#lib.reasdd(s.encode('ascii'),len(s))
print (type(b) ) 
