import IsoSpecPy
from tqdm import tqdm

t = 0.0
for x in tqdm(xrange(100000)):
    i = IsoSpecPy.Iso("C100H100N100O100")
    t += i.getTheoreticalAverageMass()
print t
