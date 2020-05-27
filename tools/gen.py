from masserstein import peptides
from string import printable
from collections import Counter

for i in range(255):
    if chr(i) in printable and not chr(i) in "\t\n" and not i in [11, 12, 13]:
        s = ["/* Code:", str(i).rjust(3), " ASCII char:", chr(i), "*/      "]
    else:
        s = ["/* Code:", str(i).rjust(3), " unprintable   */      "]

    chnosse = list("CHNOS") + ["Se"]
    s.append(", ".join(str(peptides.daminoacids.get(chr(i).upper(), Counter())[el]) for el in chnosse)+',')
    print(' '.join(s))

    
