

int_fact = 1000000000.0

def integerize(iso):
    for mass, prob in zip(iso.masses, iso.probs):
        yield (int(mass*int_fact), int(prob*int_fact))


flows_loud = True
