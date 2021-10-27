import re, sys, subprocess
results = []
search_specs = 'search_specs.txt'
multiplicity = sys.argv[1]
with open(search_specs) as specs:
    for line in specs:
        if line.find('title') >= 0:
            t = re.search('\"(.+?)\"',line)
            if t:
                name = t.group(1) + multiplicity
            break
with open(name + '/_data.txt') as data:
    for line in data:
        if line.find('DONE') < 0:
            print(line)
            results.append(line.split(';'))
results.sort(key=lambda e: e[1])
for r in results:
    r[1] = 'Energy=' + str(round(float(r[1].split('=')[1]),5))
    print(' '.join(r))
with open(name+'/_data.txt', 'w') as data:
    for r in results:
        data.write(';'.join(r))


