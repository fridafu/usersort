

a = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. God jul."
b = "The quick brown fox jumps over the lazy dog. Godt nytt år."

def lev(a,b):
	if a == "":
		return len(b)
	if a == "":
		return len(b)
	if a[-1] == b[-1]:
		cost = 0
	else:
		cost = 1
	res = min([lev(a[:-1], b)+1,
           lev(a, b[:-1])+1, 
           lev(a[:-1], b[:-1]) + cost])
    return res

stringdist = lev(a,b)
print(stringdist)