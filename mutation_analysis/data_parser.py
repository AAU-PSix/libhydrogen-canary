from random import choice
from typing import Dict, List, Tuple
from os import linesep, listdir
from os.path import isfile, join


class MutationLocationData():
    def __init__(
        self,
        id: str = None,
        runs: int = None,
        mutants_killed_count: int = None,
        mutants_survived_count: int = None,
    ) -> None:
        self.id = id
        self.runs = runs
        self.mutants_killed_count = mutants_killed_count
        self.mutants_survived_count = mutants_survived_count

    @property
    def mutant_count(self) -> int:
        return self.mutants_killed_count + self.mutants_survived_count

    @property
    def exhaustive_mutation_score(self) -> float:
        if self.mutant_count == 0: return None
        return self.mutants_killed_count / self.mutant_count

    @property
    def exhaustive_mutation_score_str(self) -> str:
        return f'{self.mutants_killed_count} / ({self.mutants_killed_count} + {self.mutants_survived_count}) = {self.exhaustive_mutation_score}'

    def sampled_mutation_score(self, weight: float) -> str:
        if self.mutant_count == 0: return None
        return (self.mutants_killed_count * weight) / (self.mutant_count * weight)

    def sampled_mutation_score_str(self, weight: float) -> str:
        return f'({self.mutants_killed_count} * {weight}) / (({self.mutants_killed_count} + {self.mutants_survived_count}) * {weight}) = {self.sampled_mutation_score(weight)}'

class OriginalLocationData():
    def __init__(
        self,
        id: str = None,
        visit_count: str = None,
    ) -> None:
        self.id = id
        self.visit_count = visit_count

class OriginalExecutionData():
    def __init__(self) -> None:
        self.locations: List[OriginalLocationData] = list()
        self.visited_locations: List[str] = list()
        self.unvisited_locations: List[str] = list()
        self.visited_nodes_count: int = None
        self.unvisited_nodes_count: int = None
        self.visited_candidates_count: int = None
        self.unvisited_candidates_count: int = None
        self.tested_mutations_count: int = None
        self.untested_mutations_count: int = None

    @property
    def location_visition_ratio(self) -> Dict[OriginalLocationData, float]:
        result: Dict[OriginalLocationData, float] = dict()
        for location in self.locations:
            if location.visit_count == 0: continue
            result[location] = location.visit_count / self.total_visitations
        return result

    @property
    def total_visitations(self) -> int:
        counter = 0
        for location in self.locations:
            counter += location.visit_count
        return counter

    @property
    def location_visited_count(self) -> int:
        return len(self.visited_locations)

    @property
    def location_unvisited_count(self) -> int:
        return len(self.unvisited_locations)

    @property
    def location_count(self) -> int:
        return self.location_visited_count + self.location_unvisited_count

    @property
    def location_coverage(self) -> float:
        return self.location_visited_count / self.location_count

    @property
    def location_coverage_str(self) -> str:
        return f'{self.location_visited_count} / ({self.location_visited_count} + {self.location_unvisited_count}) = {self.location_coverage}'

    @property
    def candidates_count(self) -> int:
        return self.visited_candidates_count + self.unvisited_candidates_count

    @property
    def candidates_coverage(self) -> float:
        return self.visited_candidates_count / self.candidates_count

    @property
    def candidates_coverate_str(self) -> float:
        return f'{self.visited_candidates_count} / ({self.visited_candidates_count} + {self.unvisited_candidates_count}) = {self.candidates_coverage}'

    @property
    def mutations_count(self) -> int:
        return self.tested_mutations_count + self.untested_mutations_count

    @property
    def mutations_coverage(self) -> float:
        return self.tested_mutations_count / self.mutations_count

    @property
    def mutations_coverage_str(self) -> float:
        return f'{self.tested_mutations_count} / ({self.tested_mutations_count} + {self.untested_mutations_count}) = {self.mutations_coverage}'

class MutationFileData():
    def __init__(self) -> None:
        self.name: str = None
        self.mutation_time: float = None
        self.total_killed: int = None
        self.total_survived: int = None
        self.original_execution_data = OriginalExecutionData()
        self.mutation_locations_data: List[MutationLocationData] = list()

    def get_mutation_location(self, location: str) -> MutationLocationData:
        for mutation in self.mutation_locations_data:
            if mutation.id == location: return mutation

    @property
    def mutants_kill_count(self) -> int:
        counter = 0
        for location in self.mutation_locations_data:
            counter += location.mutants_killed_count
        return counter

    @property
    def mutants_survival_count(self) -> int:
        counter = 0
        for location in self.mutation_locations_data:
            counter += location.mutants_survived_count
        return counter

    @property
    def mutant_count(self) -> int:
        return self.mutants_kill_count + self.mutants_survival_count

    @property
    def path_based_mutation_score(self) -> float:
        return self.mutants_kill_count / self.mutant_count

    @property
    def path_based_mutation_score_str(self) -> str:
        return f'{self.mutants_kill_count} / ({self.mutants_kill_count} + {self.mutants_survival_count}) = {self.path_based_mutation_score}'

    @property
    def exhaustive_mutation_score(self) -> float:
        return self.mutants_kill_count / (self.mutant_count + self.original_execution_data.untested_mutations_count)

    @property
    def uniform_weight(self) -> float:
        return 1 / len(self.mutation_locations_data)

    @property
    def uniform_weights(self) -> Dict[str, float]:
        result: Dict[str, float] = dict()
        for location in self.mutation_locations_data:
            result[location.id] = 1 / len(self.mutation_locations_data)
        return result

    def uniform_sample_size(self, weight: float) -> float:
        return self.mutant_count * weight

    @property
    def location_weights(self) -> Dict[str, float]:
        result: Dict[str, float] = dict()
        total_inverse = 0
        for location in self.original_execution_data.location_visition_ratio:
            ratio = self.original_execution_data.location_visition_ratio[location]
            if ratio == 0: continue
            total_inverse += 1 - ratio
        for location in self.original_execution_data.location_visition_ratio:
            ratio = self.original_execution_data.location_visition_ratio[location]
            if ratio == 0: continue
            if total_inverse == 0: result[location.id] = 0
            else: result[location.id] = (1 - ratio) / total_inverse
        return result

    def weight_sample_size(self, uniform_distribution: float, weights: Dict[str, float]) -> Dict[str, float]:
        result: Dict[str, float] = dict()
        uniform_sample_size = self.uniform_sample_size(uniform_distribution)
        for location in self.mutation_locations_data:
            if location.id not in weights: continue
            result[location.id] = min(uniform_sample_size * weights[location.id], location.mutant_count) / location.mutant_count
        return result

    @property
    def path_weighted_mutation_score(self) -> float:
        return self.weight_mutation_score(
            self.weight_sample_size(self.uniform_weight, self.location_weights)
        )
        
    @property
    def uniformly_weighted_mutation_score(self) -> float:
        return self.weight_mutation_score(self.uniform_weights)

    def weight_mutation_score(self, weights: Dict[str, float]) -> float:
        survived_mutants = 0
        killed_mutants = 0
        for location in self.mutation_locations_data:
            if location.id not in weights: continue
            weight = weights[location.id]
            survived_mutants += location.mutants_survived_count * weight
            killed_mutants += location.mutants_killed_count * weight
        if killed_mutants + survived_mutants == 0: return 1
        return killed_mutants / (killed_mutants + survived_mutants)

def parse_file(path: str) -> MutationFileData:
    data = MutationFileData()

    with open(path) as file:
        lines = file.read().split(linesep)
        idx = 0
        line = lines[idx]
        split = line.split()

        def next() -> bool:
            nonlocal idx
            nonlocal line
            nonlocal split
            nonlocal lines
            idx += 1
            if idx < len(lines):
                line = lines[idx]
                split = line.split()
                return True
            return False

        while idx < len(lines):
            if line.startswith("Mutation took") and \
                split[-1] == "seconds":
                data.mutation_time = float(split[2])
                next()
            elif "Visited count" in line:
                next() # Skip header
                while len(split) == 5:
                    data.original_execution_data.locations.append(
                        OriginalLocationData(
                            split[0], int(split[3])
                        )
                    )
                    next()
            elif "Visited locations" in line:
                for location in split[2:]:
                    data.original_execution_data.visited_locations.append(
                        location.removesuffix(",")
                    )
                next()
            elif "Unvisited locations" in line:
                for location in split[2:]:
                    data.original_execution_data.unvisited_locations.append(
                        location.removesuffix(",")
                    )
                next()
            elif "Visited nodes" in line:
                data.original_execution_data.visited_nodes_count = int(split[2])
                next()
            elif "Unvisited nodes" in line:
                data.original_execution_data.unvisited_nodes_count = int(split[2])
                next()
            elif "Visited candidates" in line:
                data.original_execution_data.visited_candidates_count = int(split[2])
                next()
            elif "Unvisited candidates" in line:
                data.original_execution_data.unvisited_candidates_count = int(split[2])
                next()
            elif "Tested mutations" in line:
                data.original_execution_data.tested_mutations_count = int(split[2])
                next()
            elif "Untested mutations" in line:
                data.original_execution_data.untested_mutations_count = int(split[2])
                next()
            elif "Location" in line and len(split) == 2:
                mutation_location_data = MutationLocationData(
                    split[1]
                )
                next()

                while "Location Killed" not in line:
                    next()

                mutation_location_data.mutants_killed_count = int(split[2])
                mutation_location_data.mutants_survived_count = int(split[4])

                if mutation_location_data.exhaustive_mutation_score is not None:
                    data.mutation_locations_data.append(
                        mutation_location_data
                    )
            else:
                next()

    return data

counter = 0

parsed_data: List[MutationFileData] = list()

files = [f for f in listdir("./") if isfile(join("./", f))]
for file in files:
    if not file.endswith(".txt"):
        continue
    mutation_name = "_".join(
        file.removesuffix(".txt") \
            .split("_")[0:-1]
    )

    data = parse_file(f'./{file}')
    data.name = mutation_name
    if data.mutant_count == 0:
        continue
    parsed_data.append(data)

parsed_data.sort(reverse=True, key=lambda x: x.original_execution_data.untested_mutations_count)

data = parsed_data[0]

#print(data.name)
#print(f"  Uniform mutation score '{data.uniformly_weighted_mutation_score}'")
#print(f"  Weighted mutation score '{data.path_weighted_mutation_score}'")

# total_inverse_ratio = 0
# for location_ratio in data.original_execution_data.location_visition_ratio:
#     location = location_ratio[0]
#     ratio = location_ratio[1]
# 
#     found = False
#     for mutation_location in data.mutation_locations_data:
#         if mutation_location.id == location.id:
#             found = True
#             break
#     if not found: continue
# 
#     inverse_ratio = 1 - ratio
#     total_inverse_ratio += inverse_ratio

# Path-based: 114,   126,   29    =  269
# 33% sample: 37,62  41,58  9,57  =  88,77
#   weighted: 22,8   25,2   17,4  =  65,4

# print("Normalized inverse ratio")
# weights: Dict[str, float] = dict()
# for location_ratio in data.original_execution_data.location_visition_ratio:
#     location = location_ratio[0]
#     ratio = location_ratio[1]
#     
#     found = False
#     for mutation_location in data.mutation_locations_data:
#         if mutation_location.id == location.id:
#             found = True
#             break
#     if not found: continue
# 
#     inverse_ratio = 1 - ratio
#     normalied_inverse_ratio = inverse_ratio / total_inverse_ratio
#     # print(f"  {location_ratio[0].id} -> {ratio}, {inverse_ratio} -> {normalied_inverse_ratio}")
#     weights[location.id] = normalied_inverse_ratio

# print(f"  Uniform sample size 20% '{data.uniform_sample_size(0.2)}'")
# adjusted_total = 0
# adjusted_weights = data.weight_sample_size(0.2, weights)
# for location in adjusted_weights:
#     adjusted_weight = adjusted_weights[location]
#     weight_mutant_count = data.get_mutation_location(location).mutant_count * adjusted_weight
#     adjusted_total += weight_mutant_count
#     print(f"    {location} -> {adjusted_weight} = {weight_mutant_count}")
# print(f"  Weight sample size '{adjusted_total}'")


# print("Mutations on locations")
# print(f"  {data.mutation_locations_data[0].id} -> {data.mutation_locations_data[0].sampled_mutation_score_str(1 / 3)}")
# print(f"  {data.mutation_locations_data[1].id} -> {data.mutation_locations_data[1].sampled_mutation_score_str(1 / 3)}")
# print(f"  {data.mutation_locations_data[2].id} -> {data.mutation_locations_data[2].sampled_mutation_score_str(1 / 3)}")
# print(f"  {data.mutation_locations_data[0].id} -> {data.mutation_locations_data[0].sampled_mutation_score_str(0.20)}")
# print(f"  {data.mutation_locations_data[1].id} -> {data.mutation_locations_data[1].sampled_mutation_score_str(0.20)}")
# print(f"  {data.mutation_locations_data[2].id} -> {data.mutation_locations_data[2].sampled_mutation_score_str(0.60)}")

# for mutation_location in data.mutation_locations_data:
#     if mutation_location.id in weights:
#         weight = weights[mutation_location.id]
#         print(f"  {mutation_location.id} -> {mutation_location.sampled_mutation_score_str(weight)}")

if False:
    points = list()
    for data in parsed_data:
        exhaustive_ms = data.exhaustive_mutation_score
        for location in data.mutation_locations_data:
            points.append("%.2f" % round(exhaustive_ms - location.exhaustive_mutation_score, 2))
        # points.append("%.2f" % round(data.path_based_mutation_score, 2))
    print(", ".join(points))

if True:
    for data in parsed_data:
        counter += 1

        line = ""
        line += f"{data.name}"
        line += f" & %.2f" % round(data.mutation_time, 2)
        line += f" & {data.original_execution_data.location_count}"
        line += f" & {data.original_execution_data.location_unvisited_count}"
        line += f" & %.2f" % round(data.original_execution_data.location_coverage, 2)
        line += f" & {data.original_execution_data.candidates_count}"
        line += f" & {data.original_execution_data.unvisited_candidates_count}"
        line += f" & %.2f" % round(data.original_execution_data.candidates_coverage, 2)
        line += f" & {data.original_execution_data.mutations_count}"
        line += f" & {data.original_execution_data.untested_mutations_count}"
        line += f" & %.2f" % round(data.original_execution_data.mutations_coverage, 2)
        
        line += f" & %.2f" % round(data.path_based_mutation_score, 2)
        line += f" & %.2f" % round(data.exhaustive_mutation_score, 2)
        # line += f" & %.2f" % round(data.path_weighted_mutation_score, 2)
        # line += f" & %.2f" % round(data.path_based_mutation_score - data.random_mutation_score, 2)
        line = line.replace("_", "\\_")
        print(line + " \\\\")

        # print()
        # print()

        # print(f"Original execution of '{mutation_name}'")
        # print(f"  Covered locations: '" + ", ".join(data.original_execution_data.visited_locations) + "'")
        # print(f"  Uncovered locations: '" + ", ".join(data.original_execution_data.unvisited_locations) + "'")
        # print(f"  total visitations: '{data.original_execution_data.total_visitations}'")
        # print(f"    Location coverage '{data.original_execution_data.location_coverage_str}'")
        # for location in data.original_execution_data.locations:
        #     print(f"      '{location.id}' was visited '{location.visit_count}' times")
        # print(f"  Candiate coverage '{data.original_execution_data.candidates_coverate_str}'")
        # print(f"  Mutation coverage '{data.original_execution_data.mutations_coverate_str}'")
        # print(f"Total mutation time {data.mutation_time} seconds")
        # print(f"  Mutation score for program was '{data.mutation_score_str}'")
        # for location in data.mutation_locations_data:
        #     print(f"Mutation at location '{location.id}'")
        #     print(f"  Mutation score for location was '{location.mutation_score_str}'")

    total_mutation_time = 0
    total_visited_locations = 0
    total_unvisited_locations = 0
    total_location_coverage = 0
    total_visited_candidate = 0
    total_unvisited_candidate = 0
    total_candidate_coverage = 0
    total_tested_mutants = 0
    total_untested_mutants = 0
    total_mutant_coverage = 0
    total_ms_path_based = 0
    total_ms_random = 0
    total_path_weighted_mutation_score = 0
    for data in parsed_data:
        total_mutation_time += data.mutation_time
        total_visited_locations += len(data.original_execution_data.visited_locations)
        total_unvisited_locations += len(data.original_execution_data.unvisited_locations)
        total_location_coverage += data.original_execution_data.location_coverage
        total_visited_candidate += data.original_execution_data.visited_candidates_count
        total_unvisited_candidate += data.original_execution_data.unvisited_candidates_count
        total_candidate_coverage += data.original_execution_data.candidates_coverage
        total_tested_mutants += data.original_execution_data.tested_mutations_count
        total_untested_mutants += data.original_execution_data.untested_mutations_count
        total_mutant_coverage += data.original_execution_data.mutations_coverage
        total_ms_path_based += data.path_based_mutation_score
        total_ms_random += data.exhaustive_mutation_score
        total_path_weighted_mutation_score += data.path_weighted_mutation_score


    # total_line = "Total"
    # total_line += f' & {total_visited_locations}'
    # total_line += f' & {total_unvisited_locations}'
    # total_line += f' & '
    # total_line += f' & {total_visited_candidate}'
    # total_line += f' & {total_unvisited_candidate}'
    # total_line += f' & '
    # total_line += f' & {total_tested_mutants}'
    # total_line += f' & {total_untested_mutants}'
    # print(total_line + "\\\\")

    average_line = "Average"
    total_amount_of_data = len(parsed_data)
    average_line += f' & %.2f ' % round(total_mutation_time / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_visited_locations / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_unvisited_locations / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_location_coverage / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_visited_candidate / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_unvisited_candidate / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_candidate_coverage / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_tested_mutants / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_untested_mutants / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_mutant_coverage / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_ms_path_based / total_amount_of_data, 2)
    average_line += f' & %.2f ' % round(total_ms_random / total_amount_of_data, 2)
    # average_line += f' & %.2f ' % round(total_path_weighted_mutation_score / total_amount_of_data, 2)
    print(average_line + "\\\\")


    print()
