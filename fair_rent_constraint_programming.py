from ortools.linear_solver import pywraplp
from itertools import permutations

solver = pywraplp.Solver('SolveIntegerProblem', pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

# Decision Variables
DV = {}

# The rent that each housemate will pay for the room assigned to them.
for h in range(1, 4):
    name = 'H'+str(h)
    DV[name] = solver.NumVar(0.0, solver.infinity(), name)

# Assigned Rooms
# These values will be 0 or 1 depending on which room a housemate is assigned to. For example, if
# housemate 1 is assigned to room 1, then H1_R1 will equal 1 and H1_R2 and H1_R3 will equal 0.
for h in range(1, 4):
    for r in range(1, 4):
        name = 'H' + str(h) + '_R' + str(r)
        DV[name] = solver.IntVar(0.0, 1.0, name)
        
# Define house mate utilities for convenience
for i in range(1, 4):
    for j in range(1, 4):
        name = 'H' + str(i) + '_swap_H' + str(j) + '_utility'
        DV[name] = solver.NumVar(-solver.infinity(), solver.infinity(), name)

# Minimum Utility, the minimum of the three house mate utilities.
DV['min_stay_utility'] = solver.NumVar(-solver.infinity(), solver.infinity(), 'min__stay_utility')

H1_preferences = [35, 50, 65]
H2_preferences = [45, 52, 53]
H3_preferences = [30, 45, 75]

# check prefences all sum to the same amount
assert sum(H1_preferences) == sum(H2_preferences) == sum(H3_preferences)

preferences = {}
H_list = [H1_preferences, H2_preferences, H3_preferences]

for index, h in enumerate(H_list):
    for r in range(0, 3):
        name = 'H' + str(index+1) + '_R' + str(r+1)
        preferences[name] = h[r]

# constraints

# The rents will sum to the total rent
# H1 + H2 + H3 == 150
c0 = solver.Constraint(sum(H1_preferences), sum(H1_preferences))
for h in range(1, 4):
    name = 'H' + str(h)
    c0.SetCoefficient(DV[name], 1)

# One room per housemate
# H1_R1 + H1_R2 + H1_R3 == 1
# H2_R1 + H2_R2 + H2_R3 == 1
# H3_R1 + H3_R2 + H3_R3 == 1
for h in range(1, 4):
    c = solver.Constraint(1, 1)
    for r in range(1, 4):
        c.SetCoefficient(DV['H' + str(h) + '_R' + str(r)], 1)

# One housemate per room
# H1_R1 + H2_R1 + H3_R1 == 1
# H1_R2 + H2_R2 + H3_R2 == 1
# H1_R3 + H2_R3 + H3_R3 == 1
for r in range(1, 4):
    c = solver.Constraint(1, 1)
    for h in range(1, 4):
        c.SetCoefficient(DV['H' + str(h) + '_R' + str(r)], 1)

# Utility
# H1_utility = H1_R1*35 + H1_R2*50 + H1_R3*65 – H1
# H2_utility = H2_R1*45 + H2_R2*52 + H2_R3*53 – H2
# H3_utility = H3_R1*30 + H3_R2*45 + H3_R3*75 – H3
# Utility of swapped rooms.
# H1_swap_H2_utility = H2_R1*35 + H2_R2*50 + H2_R3*65 – H2
# H1_swap_H3_utility = H3_R1*35 + H3_R2*50 + H3_R3*65 – H3
# H2_swap_H1_utility = H1_R1*45 + H1_R2*52 + H1_R3*53 – H1
# H2_swap_H3_utility = H3_R1*45 + H3_R2*52 + H3_R3*53 – H3
# H3_swap_H1_utility = H1_R1*30 + H1_R2*45 + H1_R3*75 – H1
# H3_swap_H2_utility = H2_R1*30 + H2_R2*45 + H2_R3*75 – H2
for h_i in range(1, 4):
    for h_j in range(1, 4):
        c = solver.Constraint(0.0, 0.0)
        name = 'H' + str(h_i) + '_swap_H' + str(h_j) + '_utility'
        c.SetCoefficient(DV[name], -1)
        c.SetCoefficient(DV['H' + str(h_j)], -1)
        for r in range(1, 4):
            c.SetCoefficient(DV['H' + str(h_j) + '_R' + str(r)], preferences['H' + str(h_i) + '_R' + str(r)])
        
# H1_swap_H2_utility, H1_swap_H3_utility <= H1_utility
# H2_swap_H1_utility, H2_swap_H3_utility <= H2_utility
# H3_swap_H1_utility, H3_swap_H2_utility <= H3_utility
perm = permutations([1, 2, 3], 2) 
for h_i, h_j in perm:
    # print(a, ' ', b)
    c = solver.Constraint(-solver.infinity(), 0.0)
    c.SetCoefficient(DV['H' + str(h_i) + '_swap_H' + str(h_i) + '_utility'], -1)
    c.SetCoefficient(DV['H' + str(h_i) + '_swap_H' + str(h_j) + '_utility'], 1)

# Make the min_stay_utility smaller or equal to the house mate
# utilities for staying
for h in range(1, 4):
    c = solver.Constraint(-solver.infinity(), 0.0)
    c.SetCoefficient(DV['min_stay_utility'], 1)
    c.SetCoefficient(DV['H' + str(h) + '_swap_H' + str(h) + '_utility'], -1)

objective = solver.Objective()
objective.SetCoefficient(DV['min_stay_utility'], 1)

objective.SetMaximization()

"""Solve the problem and print the solution."""
result_status = solver.Solve()
# The problem has an optimal solution.
assert result_status == pywraplp.Solver.OPTIMAL

# The solution looks legit (when using solvers other than
# GLOP_LINEAR_PROGRAMMING, verifying the solution is highly recommended!).
assert solver.VerifySolution(1e-7, True)

print('Number of variables =', solver.NumVariables())
print('Number of constraints =', solver.NumConstraints())

# The objective value of the solution.
print('Optimal objective value = %d' % solver.Objective().Value())
print()
# The value of each variable in the solution.


for key, value in DV.items():
    print(key, ' = ', round(value.solution_value()))
