"""
This script checks the given marxan solution's cost and shortfall.
Use this to independently validate whether a marxan solution 
"""
import argparse
import csv
import sys

def read_pu(pu_map, pu_path):
    with open(pu_path) as file:
        dialect = csv.Sniffer().sniff(file.read(), delimiters=';, \t')
        file.seek(0)
        csvread = csv.reader(file, dialect)
        headers = next(csvread)
        print(headers)

        for line in csvread:
            pu_map[int(line[0])] = float(line[1])

def read_puvspr(puvspr_map, puvspr_path):
    spec_totals = {}
    with open(puvspr_path) as file:
        dialect = csv.Sniffer().sniff(file.read(), delimiters=';, \t')
        file.seek(0)
        csvread = csv.reader(file, dialect)
        headers = next(csvread)
        print(headers)

        # ordering should be spec, pu, amount
        for line in csvread:
            puvspr_map[(int(line[0]),int(line[1]))] = float(line[2])

            if int(line[0]) in spec_totals:
                spec_totals[int(line[0])] += float(line[2])
            else:
                spec_totals[int(line[0])] = float(line[2])

    return spec_totals

def read_solution(solution_map, solution_path):
    with open(solution_path) as file:
        dialect = csv.Sniffer().sniff(file.read(), delimiters=';, \t')
        file.seek(0)
        csvread = csv.reader(file, dialect)
        headers = next(csvread)
        print(headers)

        for line in csvread:
            solution_map[int(line[0])] = int(line[1])

def read_spec(spec_map, spec_path):
    # needs to determine if target is supplied or prop.and return true for prop, false for target
    isProp = False
    with open(spec_path) as file:
        dialect = csv.Sniffer().sniff(file.read(), delimiters=';, \t')
        file.seek(0)
        csvread = csv.reader(file, dialect)
        headers = next(csvread)
        print(headers)

        # find header with prop or target
        try:
            ind_target = headers.index("prop")
            isProp = True
        except ValueError:
            ind_target = headers.index("target")

        for line in csvread:
            spec_map[int(line[0])] = float(line[ind_target])

    return isProp

def calculate_contributions(puvspr_map, solution, contributions_map):
    for term in puvspr_map:
        if solution[int(term[1])] > 0: #if PU is switched on
            if int(term[0]) in contributions_map:
                contributions_map[int(term[0])] += float(puvspr_map[term])
            else:
                contributions_map[int(term[0])] = float(puvspr_map[term])

def calculate_cost(pu_map, solution):
    cost = 0.0

    for pu in solution:
        if solution[pu] == 1 or solution[pu] == 2:
            cost += pu_map[pu]

    return cost

def calculate_shortfall(contributions, targets):
    shortfall = 0.0

    for spec in targets:
        shortfall += max(0, targets[spec] - contributions[spec])

    return shortfall

def prop_to_target(spec_map, spec_totals):
    for spec in spec_map:
        if spec in spec_totals:
            spec_map[spec] = spec_totals[spec]*spec_map[spec]
        else:
            spec_map[spec] = 0.0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='solution', type=str, default="output/outputbest.txt", help="Path to a CSV containing the solution from marxan.")
    parser.add_argument(dest='pu', type=str, default="input/pu.dat", help="Path to the pu.dat file or equivalent.")
    parser.add_argument(dest='spec', type=str, default="input/spec.dat", help="Path to the spec.dat file or equivalent.")
    parser.add_argument(dest='puvspr', type=str, default="input/puvspr.dat", help="Path to the puvspr.dat file or equivalent.")
    args = parser.parse_args()

    # declare objects
    pu_map = {} # pu -> cost
    spec_map = {} # spec -> target
    puvspr_map = {} # tuple(pu, spec) -> amount
    solution_map = {} # pu -> status
    contributions_map = {} # spec -> amount existing in solution

    # read files & compute
    read_pu(pu_map, args.pu)
    read_solution(solution_map, args.solution)
    spec_totals = read_puvspr(puvspr_map, args.puvspr)
    calculate_contributions(puvspr_map, solution_map, contributions_map)
    is_prop = read_spec(spec_map, args.spec)

    if is_prop:
        # need to convert the props into actual target amounts
        prop_to_target(spec_map, spec_totals)

    # calculate shortfall & cost
    returned_cost = calculate_cost(pu_map, solution_map)
    returned_shortfall = calculate_shortfall(contributions_map, spec_map)

    print("Calculated cost: {0} Calculated shortfall: {1}".format(returned_cost, returned_shortfall))


