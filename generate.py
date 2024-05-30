import random
import math

class City:
    def __init__(self, id, name, x, y):
        self.id = id
        self.name = name
        self.x = x
        self.y = y
        self.neighbors = []

def euclidean_distance(city1, city2):
    return math.sqrt((city1.x - city2.x) ** 2 + (city1.y - city2.y) ** 2)

def generate_cities(num_cities, min_coord, max_coord):
    cities = []
    for i in range(num_cities):
        name = chr(65 + i)  # City names A, B, C, ...
        x = random.uniform(min_coord, max_coord)
        y = random.uniform(min_coord, max_coord)
        cities.append(City(i, name, x, y))
    return cities

def calculate_distances(cities):
    for i in range(len(cities)):
        for j in range(i + 1, len(cities)):
            distance = euclidean_distance(cities[i], cities[j])
            cities[i].neighbors.append((cities[j].id, distance))
            cities[j].neighbors.append((cities[i].id, distance))

    return cities

def save_cities_to_file(cities, filename):
    with open(filename, 'w') as f:
        for city in cities:
            f.write(f"{city.id} {city.name} {city.x} {city.y}\n")
            for neighbor_id, distance in city.neighbors:
                f.write(f"{neighbor_id} {distance}\n")

# Generate random cities
num_cities = 10
min_coord = 0
max_coord = 100
cities = generate_cities(num_cities, min_coord, max_coord)

# Calculate distances between cities
cities = calculate_distances(cities)

# Save cities and their distances to a file
save_cities_to_file(cities, "input.txt")
