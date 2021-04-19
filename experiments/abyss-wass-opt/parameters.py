

int_fact = 1000000000.0
emp_grad_dval = 0.0001

assert int_fact * emp_grad_dval > 1.0

def integerize(iso):
    for mass, prob in zip(iso.masses, iso.probs):
        yield (int(mass*int_fact), int(prob*int_fact))


flows_loud = True
