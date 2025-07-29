#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <random>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <unordered_set>
#include <queue>
#include <map>

using namespace std;

// Structure to represent a customer with their details
struct Customer {
    int id;
    double x, y;
    int demand;
    double serviceTime;
    double earliest;
    double latest;
};

// Structure to represent a route with customer IDs, distances, and times
struct Route {
    vector<int> customerIds;
    double totalDistance;
    vector<double> arrivalTimes;
    vector<double> departureTimes;
    int load; // Total demand of customers in the route
    Route() : totalDistance(0.0), load(0) {}
};

// Structure of problem data
struct ProblemData {
    vector<Customer> customers;
    int vehicleCapacity;
    int maxVehicles;
    double maxRouteDuration;
};

// Euclidean distance between two customers
double distance(const Customer &a, const Customer &b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

// Reads problem instance
ProblemData readInstance(const string &filename) {
    ProblemData data;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Cannot open file: " << filename << endl;
        exit(1);
    }
    string line;

    // Skip header lines and read vehicle info
    while (getline(infile, line)) {
        if (line.find("CUST NO.") != string::npos) {
            getline(infile, line); // Skip the empty line
            break;
        } else if (line.find("NUMBER") != string::npos) {
            getline(infile, line);
            istringstream iss(line);
            iss >> data.maxVehicles >> data.vehicleCapacity;
        }
    }

    // Read customer data
    int numCustomers = 0;
    while (getline(infile, line)) {
        istringstream issCust(line);
        Customer cust;
        if (issCust >> cust.id >> cust.x >> cust.y >> cust.demand >> cust.earliest >> cust.latest >> cust.serviceTime) {
            data.customers.push_back(cust);
            numCustomers++;
        } else {
            break; // Stop if line is invalid
        }
    }
    data.maxRouteDuration = data.customers[0].latest; // Set from depot's due time
    infile.close();
    return data;
}

// Checks if a single route is feasible based on capacity, time windows, and duration
bool isRouteFeasible(const Route &route, const ProblemData &data) {
    int capacityUsed = 0;
    double currentTime = 0.0;
    int currentIndex = 0; // Start from depot (index 0)

    for (int custId : route.customerIds) {
        // Find customer index by ID
        auto it = find_if(data.customers.begin(), data.customers.end(), [custId](const Customer &c) {
            return c.id == custId;
        });
        if (it == data.customers.end()) {
            cerr << "Customer with id " << custId << " not found in isRouteFeasible." << endl;
            return false;
        }
        int index = distance(data.customers.begin(), it);

        // Capacity check
        capacityUsed += data.customers[index].demand;
        if (capacityUsed > data.vehicleCapacity) {
            return false;
        }

        // Time window check
        double travelTime = distance(data.customers[currentIndex], data.customers[index]);
        double arrivalTime = currentTime + travelTime;
        double serviceStartTime = max(arrivalTime, data.customers[index].earliest);
        if (serviceStartTime > data.customers[index].latest) {
            return false;
        }
        currentTime = serviceStartTime + data.customers[index].serviceTime;
        currentIndex = index; // Update current index for next customer
    }

    // Return to depot
    double returnTravelTime = distance(data.customers[currentIndex], data.customers[0]);
    currentTime += returnTravelTime;
    if (currentTime > data.maxRouteDuration) {
        return false;
    }
    return true;
}

// Checks if the entire solution is feasible
bool isSolutionFeasible(const vector<Route>& routes, const ProblemData& data) {
    // 1. Check feasibility of each route
    for (const auto& route : routes) {
        if (!isRouteFeasible(route, data)) {
            return false;
        }
    }

    // 2. Check if all customers are visited exactly once
    vector<bool> visited(data.customers.size(), false);
    visited[0] = true; // Mark depot as visited
    for (const auto& route : routes) {
        for (int custId : route.customerIds) {
            auto it = find_if(data.customers.begin(), data.customers.end(), [custId](const Customer &c) {
                return c.id == custId;
            });
            int index = distance(data.customers.begin(), it);
            if (visited[index]) {
                return false; // Customer visited more than once
            }
            visited[index] = true;
        }
    }

    // Check if all customers (except depot) are visited
    for (size_t i = 1; i < visited.size(); ++i) {
        if (!visited[i]) {
            return false;
        }
    }
    return true;
}

// Finds the index of a customer in the customers list by their ID
int find_customer_index(const vector<Customer>& customers, int custId) {
    for (size_t i = 0; i < customers.size(); ++i) {
        if (customers[i].id == custId) {
            return i;
        }
    }
    return -1; // Not found
}

// Constructs an initial solution using a greedy nearest neighbor approach
vector<Route> constructInitialSolution(const ProblemData &data) {
    vector<Route> routes;
    vector<bool> visited(data.customers.size(), false);
    visited[0] = true; // Depot is always visited

    while (find(visited.begin() + 1, visited.end(), false) != visited.end()) {
        Route route;
        route.arrivalTimes.push_back(0.0); // Depot arrival time (always 0)
        double currentTime = 0.0;
        int currentId = 0; // Start from the depot
        int currentLoad = 0;

        while (true) {
            int nextCustomer = -1; // Initialize to -1
            double bestDist = numeric_limits<double>::max(); // Initialize to infinity

            for (size_t i = 1; i < data.customers.size(); i++) { // Start from 1 (skip depot)
                if (!visited[i]) {
                    double travelTime = distance(data.customers[currentId], data.customers[i]);
                    double arrivalTime = currentTime + travelTime;
                    if (arrivalTime < data.customers[i].earliest) {
                        arrivalTime = data.customers[i].earliest;
                    }
                    if (arrivalTime <= data.customers[i].latest) {
                        if (travelTime < bestDist) {
                            bestDist = travelTime;
                            nextCustomer = i;
                        }
                    }
                }
            }

            if (nextCustomer == -1) break; // No more feasible customers

            // Calculate arrival and departure times
            double arrivalTime = currentTime + distance(data.customers[currentId], data.customers[nextCustomer]);
            if (arrivalTime < data.customers[nextCustomer].earliest) {
                arrivalTime = data.customers[nextCustomer].earliest;
            }
            double departureTime = arrivalTime + data.customers[nextCustomer].serviceTime;
            int newLoad = currentLoad + data.customers[nextCustomer].demand;

            // Check capacity
            if (newLoad <= data.vehicleCapacity) {
                Route tempRoute = route;
                tempRoute.customerIds.push_back(data.customers[nextCustomer].id);
                tempRoute.arrivalTimes.push_back(arrivalTime);
                tempRoute.departureTimes.push_back(departureTime);
                tempRoute.load = newLoad;

                // Check full feasibility (time windows)
                if (isRouteFeasible(tempRoute, data)) {
                    route = tempRoute; // Copy the temporary route to the actual route
                    currentTime = departureTime; // Update current time
                    visited[nextCustomer] = true; // Mark as visited
                    currentId = data.customers[nextCustomer].id;
                    currentLoad = newLoad; // Update load
                } else {
                    break; // This customer is infeasible, try a new route
                }
            } else {
                break; // Capacity exceeded, try a new route
            }
        }

        // Add return-to-depot information only if the route is not empty
        if (!route.customerIds.empty()) {
            int lastCustomerIndex = find_customer_index(data.customers, currentId);
            double returnTime = currentTime + distance(data.customers[lastCustomerIndex], data.customers[0]);
            route.arrivalTimes.push_back(returnTime);
            route.departureTimes.push_back(returnTime);
            route.totalDistance += distance(data.customers[lastCustomerIndex], data.customers[0]);
            routes.push_back(route);
        }
    }
    return routes;
}

//Random
// vector<Route> constructInitialSolution(const ProblemData &data) {
//     vector<Route> routes;
//     vector<bool> visited(data.customers.size(), false);
//     visited[0] = true; // Depot is always visited

//     //initialize random number generator
//     std::mt19937 rng(static_cast<unsigned>(std::time(nullptr)));

//     // Continue until all customers (except depot) are visited
//     while (find(visited.begin() + 1, visited.end(), false) != visited.end()) {
//         Route route;
//         route.arrivalTimes.push_back(0.0); // Depot arrival time (always 0)
//         double currentTime = 0.0;
//         int currentId = 0; // Start from the depot
//         int currentLoad = 0;

//         while (true) {
//             // Create a list of feasible customers (not visited and satisfying constraints)
//             vector<int> feasibleCustomers;
//             for (size_t i = 1; i < data.customers.size(); i++) { // Start from 1 (skip depot)
//                 if (!visited[i]) {
//                     double travelTime = distance(data.customers[currentId], data.customers[i]);
//                     double arrivalTime = currentTime + travelTime;
//                     // If arriving too early, wait until the earliest time
//                     if (arrivalTime < data.customers[i].earliest) {
//                         arrivalTime = data.customers[i].earliest;
//                     }
//                     // Check time window constraint
//                     if (arrivalTime <= data.customers[i].latest) {
//                         int newLoad = currentLoad + data.customers[i].demand;
//                         // Check capacity constraint
//                         if (newLoad <= data.vehicleCapacity) {
//                             Route tempRoute = route;
//                             tempRoute.customerIds.push_back(data.customers[i].id);
//                             tempRoute.arrivalTimes.push_back(arrivalTime);
//                             tempRoute.departureTimes.push_back(arrivalTime + data.customers[i].serviceTime);
//                             tempRoute.load = newLoad;
//                             // Check full feasibility (time windows and other constraints)
//                             if (isRouteFeasible(tempRoute, data)) {
//                                 feasibleCustomers.push_back(i); // Add to feasible customers
//                             }
//                         }
//                     }
//                 }
//             }

//             // If no feasible customers are found, end this route
//             if (feasibleCustomers.empty()) break;

//             // Randomly select a customer from the feasible ones
//             std::uniform_int_distribution<size_t> customerDist(0, feasibleCustomers.size() - 1);
//             int nextCustomer = feasibleCustomers[customerDist(rng)];

//             // Calculate arrival and departure times *only once*
//             double arrivalTime = currentTime + distance(data.customers[currentId], data.customers[nextCustomer]);
//             if (arrivalTime < data.customers[nextCustomer].earliest) {
//                 arrivalTime = data.customers[nextCustomer].earliest;
//             }
//             double departureTime = arrivalTime + data.customers[nextCustomer].serviceTime;
//             int newLoad = currentLoad + data.customers[nextCustomer].demand;

//             // Add the selected customer to the route
//             route.customerIds.push_back(data.customers[nextCustomer].id);
//             route.arrivalTimes.push_back(arrivalTime);
//             route.departureTimes.push_back(departureTime);
//             route.load = newLoad;

//             // Update current state
//             currentTime = departureTime;
//             visited[nextCustomer] = true;
//             currentId = data.customers[nextCustomer].id;
//             currentLoad = newLoad;
//         } // End of inner while loop (building a single route)

//         // Add return-to-depot information *only if the route is not empty*
//         if (!route.customerIds.empty()) {
//             int lastCustomerIndex = find_customer_index(data.customers, currentId);
//             double returnTime = currentTime + distance(data.customers[lastCustomerIndex], data.customers[0]);
//             route.arrivalTimes.push_back(returnTime);
//             route.departureTimes.push_back(returnTime);
//             route.totalDistance += distance(data.customers[lastCustomerIndex], data.customers[0]);
//             routes.push_back(route);
//         }
//     } // End of outer while loop (building all routes)

//     return routes;
// }


// Objective function: number of vehicles and total distance
pair<double, double> objectiveFunction(const vector<Route> &routes) {
    int vehicles = 0;
    double totalDistance = 0.0;
    for (const auto &route : routes) {
        if (!route.customerIds.empty()) {
            vehicles++;
            totalDistance += route.totalDistance;
        }
    }
    return {vehicles, totalDistance};
}

// Compares two solutions to determine if the new one is better
bool BetterSolution(const vector<Route> &newRoutes, const vector<Route> &currentRoutes) {
    auto [vehiclesNew, distNew] = objectiveFunction(newRoutes);
    auto [vehiclesCurrent, distCurrent] = objectiveFunction(currentRoutes);
    if (vehiclesNew < vehiclesCurrent) return true;
    if (vehiclesNew > vehiclesCurrent) return false;
    return distNew < distCurrent;
}

// Function to update route times and load after a change
void updateRoute(Route& route, const ProblemData& data) {
    route.arrivalTimes.clear();
    route.departureTimes.clear();
    route.load = 0;
    route.totalDistance = 0;
    if (route.customerIds.empty()) {
        return;
    }
    double currentTime = 0.0;
    route.arrivalTimes.push_back(currentTime); // Start time with depot
    int currentIndex = 0;
    for (int custId : route.customerIds) {
        if (custId < 0 || custId >= static_cast<int>(data.customers.size())) {
            cerr << "Invalid customer ID " << custId << " in updateRoute!" << endl;
            return;
        }
        double travelTime = distance(data.customers[currentIndex], data.customers[custId]);
        currentTime += travelTime;
        route.totalDistance += travelTime;
        if (currentTime < data.customers[custId].earliest) {
            currentTime = data.customers[custId].earliest;
        }
        route.arrivalTimes.push_back(currentTime);
        currentTime += data.customers[custId].serviceTime;
        route.departureTimes.push_back(currentTime);
        route.load += data.customers[custId].demand;
        currentIndex = custId;
    }
    // Return to depot
    double return_dist = distance(data.customers[currentIndex], data.customers[0]);
    route.totalDistance += return_dist;
    currentTime += return_dist;
    route.arrivalTimes.push_back(currentTime);
    route.departureTimes.push_back(currentTime);
}

// Validates if a customer ID is within the valid range
bool isValidCustomerId(int custId, const ProblemData& data) {
    return custId >= 0 && custId < data.customers.size();
}

// Removes a route and redistributes its customers to other routes
void removeRouteMove(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    if (routes.size() <= 1) {
        return;
    }
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t worstRouteIdx = routeDist(rng);
    std::vector<int> customersToRedistribute = routes[worstRouteIdx].customerIds;
    routes.erase(routes.begin() + worstRouteIdx);

    for (int customer : customersToRedistribute) {
        bool inserted = false;
        std::vector<size_t> routeIndices(routes.size());
        std::iota(routeIndices.begin(), routeIndices.end(), 0);
        std::shuffle(routeIndices.begin(), routeIndices.end(), rng);

        for (size_t r : routeIndices) {
            if (inserted) break;
            std::vector<size_t> positions(routes[r].customerIds.size() + 1);
            std::iota(positions.begin(), positions.end(), 0);
            std::shuffle(positions.begin(), positions.end(), rng);
            for (size_t pos : positions) {
                routes[r].customerIds.insert(routes[r].customerIds.begin() + pos, customer);
                updateRoute(routes[r], data);
                if (isRouteFeasible(routes[r], data)) {
                    inserted = true;
                    break;
                } else {
                    routes[r].customerIds.erase(routes[r].customerIds.begin() + pos);
                    updateRoute(routes[r], data);
                }
            }
        }

        if (!inserted) {
            Route newRoute;
            newRoute.customerIds.push_back(customer);
            updateRoute(newRoute, data);
            if (isRouteFeasible(newRoute, data)) {
                routes.push_back(newRoute);
            }
        }
    }

    routes.erase(std::remove_if(routes.begin(), routes.end(),
                                [](const Route& r) { return r.customerIds.empty(); }),
                 routes.end());
}

void twoOpt(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    if (routes.empty()) return;
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t routeIdx = routeDist(rng);
    while (routes[routeIdx].customerIds.empty()) {
        routeIdx = routeDist(rng);
    }
    Route& route = routes[routeIdx];
    if (route.customerIds.size() < 4) return;

    std::uniform_int_distribution<size_t> posDist(1, route.customerIds.size() - 3);
    size_t i = posDist(rng);
    size_t j = i + 1 + posDist(rng) % (route.customerIds.size() - i - 2);
    if (j >= route.customerIds.size()) j = route.customerIds.size() - 1;

    std::vector<int> originalCustomers = route.customerIds;
    std::reverse(route.customerIds.begin() + i, route.customerIds.begin() + j + 1);
    updateRoute(route, data);

    if (!isRouteFeasible(route, data)) {
        route.customerIds = originalCustomers;
        updateRoute(route, data);
    }
}

void relocateCustomer(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    if (routes.size() < 2) {
        return;
    }
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t sourceRouteIdx = routeDist(rng);
    while (routes[sourceRouteIdx].customerIds.empty()) {
        sourceRouteIdx = routeDist(rng);
    }
    std::uniform_int_distribution<size_t> custDist(0, routes[sourceRouteIdx].customerIds.size() - 1);
    size_t custPos = custDist(rng);
    int customer = routes[sourceRouteIdx].customerIds[custPos];

    routes[sourceRouteIdx].customerIds.erase(routes[sourceRouteIdx].customerIds.begin() + custPos);
    updateRoute(routes[sourceRouteIdx], data);

    size_t destRouteIdx = routeDist(rng);
    while (destRouteIdx == sourceRouteIdx) {
        destRouteIdx = routeDist(rng);
    }
    bool inserted = false;

    std::vector<size_t> positions(routes[destRouteIdx].customerIds.size() + 1);
    std::iota(positions.begin(), positions.end(), 0);
    std::shuffle(positions.begin(), positions.end(), rng);

    for (size_t pos : positions) {
        routes[destRouteIdx].customerIds.insert(routes[destRouteIdx].customerIds.begin() + pos, customer);
        updateRoute(routes[destRouteIdx], data);
        if (isRouteFeasible(routes[destRouteIdx], data)) {
            inserted = true;
            break;
        } else {
            routes[destRouteIdx].customerIds.erase(routes[destRouteIdx].customerIds.begin() + pos);
            updateRoute(routes[destRouteIdx], data);
        }
    }

    if (!inserted) {
        std::vector<size_t> routeIndices(routes.size());
        std::iota(routeIndices.begin(), routeIndices.end(), 0);
        std::shuffle(routeIndices.begin(), routeIndices.end(), rng);
        for (size_t r : routeIndices) {
            if (r == sourceRouteIdx) continue;
            if (inserted) break;
            std::vector<size_t> positions(routes[r].customerIds.size() + 1);
            std::iota(positions.begin(), positions.end(), 0);
            std::shuffle(positions.begin(), positions.end(), rng);
            for (size_t pos : positions) {
                routes[r].customerIds.insert(routes[r].customerIds.begin() + pos, customer);
                updateRoute(routes[r], data);
                if (isRouteFeasible(routes[r], data)) {
                    inserted = true;
                    break;
                } else {
                    routes[r].customerIds.erase(routes[r].customerIds.begin() + pos);
                    updateRoute(routes[r], data);
                }
            }
        }
    }

    if (!inserted) {
        Route newRoute;
        newRoute.customerIds.push_back(customer);
        updateRoute(newRoute, data);
        if (isRouteFeasible(newRoute, data)) {
            routes.push_back(newRoute);
        }
    }

    routes.erase(std::remove_if(routes.begin(), routes.end(),
                                [](const Route& r) { return r.customerIds.empty(); }),
                 routes.end());
}

void insertionMove(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    if (routes.empty()) {
        return;
    }
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t sourceRouteIdx = routeDist(rng);
    while (routes[sourceRouteIdx].customerIds.empty()) {
        sourceRouteIdx = routeDist(rng);
    }
    std::uniform_int_distribution<size_t> custDist(0, routes[sourceRouteIdx].customerIds.size() - 1);
    size_t custPos = custDist(rng);
    int customer = routes[sourceRouteIdx].customerIds[custPos];

    routes[sourceRouteIdx].customerIds.erase(routes[sourceRouteIdx].customerIds.begin() + custPos);
    updateRoute(routes[sourceRouteIdx], data);

    double minInsertionCost = std::numeric_limits<double>::max();
    size_t bestRouteIdx = -1;
    size_t bestPos = -1;

    for (size_t r = 0; r < routes.size(); ++r) {
        if (r == sourceRouteIdx) continue;
        for (size_t pos = 0; pos <= routes[r].customerIds.size(); ++pos) {
            double insertionCost = 0.0;
            if (routes[r].customerIds.empty()) {
                insertionCost = 2 * distance(data.customers[0], data.customers[customer]);
            } else {
                int prevCust = (pos == 0) ? 0 : routes[r].customerIds[pos - 1];
                int nextCust = (pos == routes[r].customerIds.size()) ? 0 : routes[r].customerIds[pos];
                insertionCost = distance(data.customers[prevCust], data.customers[customer]) +
                                distance(data.customers[customer], data.customers[nextCust]) -
                                distance(data.customers[prevCust], data.customers[nextCust]);
            }
            routes[r].customerIds.insert(routes[r].customerIds.begin() + pos, customer);
            updateRoute(routes[r], data);
            if (isRouteFeasible(routes[r], data) && insertionCost < minInsertionCost) {
                minInsertionCost = insertionCost;
                bestRouteIdx = r;
                bestPos = pos;
            }
            routes[r].customerIds.erase(routes[r].customerIds.begin() + pos);
            updateRoute(routes[r], data);
        }
    }

    bool inserted = false;
    if (bestRouteIdx != -1) {
        routes[bestRouteIdx].customerIds.insert(routes[bestRouteIdx].customerIds.begin() + bestPos, customer);
        updateRoute(routes[bestRouteIdx], data);
        inserted = true;
    }

    if (!inserted) {
        Route newRoute;
        newRoute.customerIds.push_back(customer);
        updateRoute(newRoute, data);
        if (isRouteFeasible(newRoute, data)) {
            routes.push_back(newRoute);
        }
    }

    routes.erase(std::remove_if(routes.begin(), routes.end(),
                                [](const Route& r) { return r.customerIds.empty(); }),
                 routes.end());
}

void swapCustomers(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    if (routes.size() < 2) {
        return;
    }
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t route1Idx = routeDist(rng);
    while (routes[route1Idx].customerIds.empty()) {
        route1Idx = routeDist(rng);
    }
    size_t route2Idx = routeDist(rng);
    while (route2Idx == route1Idx || routes[route2Idx].customerIds.empty()) {
        route2Idx = routeDist(rng);
    }

    std::uniform_int_distribution<size_t> custDist1(0, routes[route1Idx].customerIds.size() - 1);
    size_t pos1 = custDist1(rng);
    int customer1 = routes[route1Idx].customerIds[pos1];

    std::uniform_int_distribution<size_t> custDist2(0, routes[route2Idx].customerIds.size() - 1);
    size_t pos2 = custDist2(rng);
    int customer2 = routes[route2Idx].customerIds[pos2];

    routes[route1Idx].customerIds[pos1] = customer2;
    routes[route2Idx].customerIds[pos2] = customer1;

    updateRoute(routes[route1Idx], data);
    updateRoute(routes[route2Idx], data);

    if (!isRouteFeasible(routes[route1Idx], data) || !isRouteFeasible(routes[route2Idx], data)) {
        routes[route1Idx].customerIds[pos1] = customer1;
        routes[route2Idx].customerIds[pos2] = customer2;
        updateRoute(routes[route1Idx], data);
        updateRoute(routes[route2Idx], data);
    }
}

// CROSS Exchange neighborhood
void crossExchange(vector<Route>& routes, mt19937& rng, const ProblemData& data, int maxLength) {
    if (routes.size() < 2) return;
    uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t route1Idx = routeDist(rng);
    while (routes[route1Idx].customerIds.empty()) {
        route1Idx = routeDist(rng);
    }
    size_t route2Idx = routeDist(rng);
    while (route2Idx == route1Idx || routes[route2Idx].customerIds.empty()) {
        route2Idx = routeDist(rng);
    }

    Route& route1 = routes[route1Idx];
    Route& route2 = routes[route2Idx];

    uniform_int_distribution<size_t> lenDist(1, maxLength);
    size_t len1 = min(lenDist(rng), route1.customerIds.size());
    size_t len2 = min(lenDist(rng), route2.customerIds.size());

    uniform_int_distribution<size_t> posDist1(0, route1.customerIds.size() - len1);
    size_t pos1 = posDist1(rng);
    uniform_int_distribution<size_t> posDist2(0, route2.customerIds.size() - len2);
    size_t pos2 = posDist2(rng);

    // Store original routes
    vector<int> original1 = route1.customerIds;
    vector<int> original2 = route2.customerIds;

    // Perform CROSS exchange
    vector<int> seq1(route1.customerIds.begin() + pos1, route1.customerIds.begin() + pos1 + len1);
    vector<int> seq2(route2.customerIds.begin() + pos2, route2.customerIds.begin() + pos2 + len2);

    route1.customerIds.erase(route1.customerIds.begin() + pos1, route1.customerIds.begin() + pos1 + len1);
    route1.customerIds.insert(route1.customerIds.begin() + pos1, seq2.begin(), seq2.end());
    route2.customerIds.erase(route2.customerIds.begin() + pos2, route2.customerIds.begin() + pos2 + len2);
    route2.customerIds.insert(route2.customerIds.begin() + pos2, seq1.begin(), seq1.end());

    updateRoute(route1, data);
    updateRoute(route2, data);

    // Check feasibility
    if (!isRouteFeasible(route1, data) || !isRouteFeasible(route2, data)) {
        route1.customerIds = original1;
        route2.customerIds = original2;
        updateRoute(route1, data);
        updateRoute(route2, data);
    }
}

// Route Merge Move neighborhood
void routeMergeMove(vector<Route>& routes, mt19937& rng, const ProblemData& data) {
    if (routes.size() < 2) return; // Need at least two routes to merge

    // Select two distinct routes randomly
    uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t route1Idx = routeDist(rng);
    while (routes[route1Idx].customerIds.empty()) {
        route1Idx = routeDist(rng);
    }
    size_t route2Idx = routeDist(rng);
    while (route2Idx == route1Idx || routes[route2Idx].customerIds.empty()) {
        route2Idx = routeDist(rng);
    }

    Route& route1 = routes[route1Idx];
    Route& route2 = routes[route2Idx];

    // Store original routes for rollback if needed
    vector<int> original1 = route1.customerIds;
    vector<int> original2 = route2.customerIds;

    // Merge route2 into route1 by appending route2's customers to route1
    route1.customerIds.insert(route1.customerIds.end(), route2.customerIds.begin(), route2.customerIds.end());
    route2.customerIds.clear(); // Empty route2

    // Update route1 to compute new timings, load, and feasibility
    updateRoute(route1, data);
    updateRoute(route2, data);

    // Check feasibility of the merged route
    if (!isRouteFeasible(route1, data) || route1.load > data.vehicleCapacity) {
        // If infeasible, rollback the merge
        route1.customerIds = original1;
        route2.customerIds = original2;
        updateRoute(route1, data);
        updateRoute(route2, data);
    } else {
        // Remove empty route2 from the solution
        routes.erase(routes.begin() + route2Idx);
    }
}

size_t hashSolution(const vector<Route>& routes) {
    size_t hashValue = 0;
    for (const auto& route : routes) {
        for (int custId : route.customerIds) {
            hashValue = hashValue * 31 + custId;
        }
    }
    return hashValue;
}

vector<Route> tabuSearch(const ProblemData &data,  vector<Route> initRoutes,int maxIterations,int maxTime, int maxEvaluations) {
    vector<Route> currentRoutes = initRoutes;
    vector<Route> bestRoutes = initRoutes;
    if(!currentRoutes.empty()){
        for (auto& route : currentRoutes) updateRoute(route, data);
        for (auto& route : bestRoutes) updateRoute(route, data);}
    // ShortTerm memory:store solution hashes
    unordered_set<size_t> tabuSet; //For quick lookup of existing solutions
    queue<size_t> tabuQueue; // To manage the order and remove the oldest solution

    //Long Term memory:store the frequency of moves
    map<pair<int, int>, int> moveFrequency;

    int evaluations = 0;
    auto startTime = chrono::steady_clock::now();
    mt19937 rng(static_cast<unsigned>(time(nullptr)));
    const int printInterval = 100;
    const double improvementThreshold = 0.05;

    // fixed weights for operators
    const vector<double> operator_weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0};
    bool allZero = true;
    for (double weight : operator_weights) {
        if (weight > 0.0) {
            allZero = false;
            break;}}
    if (allZero) {
        cout << "Error: All operator weights are zero" << endl;
        cout << "Returning initial solution without search." << endl;
        auto [vehicles, distance] = objectiveFunction(bestRoutes);
        cout << "Initial Solution - Vehicles: " << vehicles << " | Distance: " << fixed << setprecision(2) << distance << endl;
        return bestRoutes;}
    // Calculate tabuTenure based on problem size
    int n = data.customers.size() - 1; // Number of customers 
    int tabuTenure = static_cast<int>(2 * sqrt(n)); 
    int noImprovementCount = 0; 
    const int noImprovementLimit = 1000000; 
    discrete_distribution<int> opDist(operator_weights.begin(), operator_weights.end());
    cout << "Fixed Weights: {";
    for (size_t i = 0; i < operator_weights.size(); ++i) {
        cout << operator_weights[i] << (i == operator_weights.size() - 1 ? "" : ", ");}
    cout << "}" << endl;
    cout << "Calculated Tabu Tenure: " << tabuTenure << endl; // Print the calculated value
    cout << "------------------------" << endl;
    while (evaluations < maxIterations && (maxEvaluations == 0 || evaluations < maxEvaluations)) {
        auto currentTime = chrono::steady_clock::now();
        double elapsedTime = chrono::duration_cast<chrono::seconds>(currentTime - startTime).count();
        if (maxTime > 0 && elapsedTime >= maxTime) break;

        // Stopping condition based on lack of improvement
        if (noImprovementCount >= noImprovementLimit) {
            cout << noImprovementCount << " NoImprovement Limit" << endl;
            break;}
        vector<Route> bestNeighbor;
        auto [bestNeighborVehicles, bestNeighborDistance] = objectiveFunction(bestRoutes);
        bool foundBetter = false;
        pair<int, int> bestMove = {-1, -1};

        // select operator based on fixed weights
        int op = opDist(rng);

        vector<Route> neighborRoutes = currentRoutes;
        pair<int, int> move = {-1, -1};
        if (op == 0) {
            removeRouteMove(neighborRoutes, rng, data);
            move = {neighborRoutes[0].customerIds.empty() ? -1 : neighborRoutes[0].customerIds[0], 0};
        } else if (op == 1) {
            relocateCustomer(neighborRoutes, rng, data);
            move = {neighborRoutes[0].customerIds.empty() ? -1 : neighborRoutes[0].customerIds[0], 0};
        } else if (op == 2) {
            insertionMove(neighborRoutes, rng, data);
            move = {neighborRoutes[0].customerIds.empty() ? -1 : neighborRoutes[0].customerIds[0], 0};
        } else if (op == 3) {
            swapCustomers(neighborRoutes, rng, data);
            move = {neighborRoutes[0].customerIds.empty() ? -1 : neighborRoutes[0].customerIds[0], 0};
        } else if (op == 4) {
            twoOpt(neighborRoutes, rng, data);
            move = {-1, -1};}
             else if (op == 5) {
            crossExchange(neighborRoutes, rng, data, 8);
            move = {neighborRoutes[0].customerIds.empty() ? -1 : neighborRoutes[0].customerIds[0], 0};}
             else {
            routeMergeMove(neighborRoutes, rng, data);
            move = {neighborRoutes[0].customerIds.empty() ? -1 : neighborRoutes[0].customerIds[0], 0};}
        evaluations++;
        if (isSolutionFeasible(neighborRoutes, data)) {
            auto [neighborVehicles, neighborDistance] = objectiveFunction(neighborRoutes);
            double neighborScore = neighborVehicles * 1000.0 + neighborDistance;
            //apply penalty from Long Term memory
            int frequency = moveFrequency[move];
            double penalty = frequency * 100.0;
            double neighborScoreWithPenalty = neighborScore + penalty;
            //check ShortTerm memory
            size_t neighborHash = hashSolution(neighborRoutes);
            bool isTabu = tabuSet.find(neighborHash) != tabuSet.end();

            auto [bestVehicles, bestDistance] = objectiveFunction(bestRoutes);
            double bestScore = bestVehicles * 1000.0 + bestDistance;
            bool significantImprovement = (bestScore - neighborScoreWithPenalty) / bestScore >= improvementThreshold;
            bool canAcceptTabu = BetterSolution(neighborRoutes, bestRoutes)||significantImprovement;

            if (!isTabu || canAcceptTabu) {
                if (BetterSolution(neighborRoutes, currentRoutes)) {
                    bestNeighbor = neighborRoutes;
                    bestNeighborVehicles = neighborVehicles;
                    bestNeighborDistance = neighborDistance;
                    bestMove = move;
                    foundBetter = true;
                    if (isTabu && canAcceptTabu) {
                        cout << "Solution was tabu but accepted. ";
                    }}}}
        if (foundBetter) {
            currentRoutes = bestNeighbor;
            if (BetterSolution(currentRoutes, bestRoutes)) {
                bestRoutes = currentRoutes;
                noImprovementCount = 0;
            } else {
                noImprovementCount++;}
            // Update ShortTerm memory
            size_t currentHash = hashSolution(currentRoutes);
            tabuSet.insert(currentHash);
            tabuQueue.push(currentHash);
            if (tabuQueue.size() > tabuTenure) {
                size_t oldestHash = tabuQueue.front();
                tabuQueue.pop();
                tabuSet.erase(oldestHash);}
            //Update Longterm memory
            if (bestMove.first != -1) {
                moveFrequency[bestMove]++;}
        } else {
            noImprovementCount++; // increment if no better solution is found
        }
        if (evaluations % printInterval == 0) {
            auto [currentBestVehicles, currentBestDistance] = objectiveFunction(bestRoutes);
            cout << "Evaluations: " << evaluations
                 << " | Vehicles: " << currentBestVehicles
                 << " | Distance: " << fixed << setprecision(0) << currentBestDistance
                 << endl;}}
    cout << "\n-----------------------------------------" << endl;
    cout << "Tabu Search Finished!" << endl;
    cout << "Total Evaluations: " << evaluations << endl;
    cout << "-----------------------------------------" << endl;
    return bestRoutes;}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << endl;
        return 1;}
    string instanceFile = argv[1];
    int maxTime = 0;
    int maxEvaluations = 0;
    try {
        maxTime = stoi(argv[2]);
        maxEvaluations = stoi(argv[3]);
        if (maxTime < 0 || maxEvaluations < 0) {
            cerr << "error:max_time/max_evaluations cannot be negative." << endl;
            return 1;}
    } catch (const std::exception& e) {
        cerr << "Error parsing args: " << e.what() << endl;
        return 1;}
    cout << "reading instance file: " << instanceFile << ".." << endl;
    ProblemData data = readInstance(instanceFile);
    if (data.customers.empty()) {
        cerr << "Error reading instance." << endl;
        return 1;}
    cout << "Instance read successfully. Customers: " << data.customers.size() - 1 << endl;
    vector<Route> initRoutes = constructInitialSolution(data);
    if (!initRoutes.empty() && !isSolutionFeasible(initRoutes, data)) {
        cerr << " Warning:Initial solution infeasible!" << endl;
    } else if (initRoutes.empty() && data.customers.size() > 1) {
        cerr << "Warning: Initial solution construction failed/empty." << endl;}

    
    cout << "\n-----------------------------------------" << endl;
    cout << "Initial Solution:" << endl;
    if (initRoutes.empty() && data.customers.size() > 1) {
        cout << "(No routes in initial solution)" << endl;} else {
        int routeNumber = 1;
        for (size_t i = 0; i < initRoutes.size(); ++i) {
            if (!initRoutes[i].customerIds.empty()) {
                cout << "Route " << routeNumber++ << ":";
                for (int custId : initRoutes[i].customerIds) {
                    cout << " " << custId;}
                cout << endl;}}
        auto [initialVehicles, initialDistance] = objectiveFunction(initRoutes);
        cout << "Vehicles: " << initialVehicles << endl;
        cout << "Distance: " << fixed << setprecision(2) << initialDistance << endl;}
    cout << "-----------------------------------------" << endl;

    int maxIterations = 10000;

    auto mainStartTime = chrono::steady_clock::now();
    vector<Route> bestRoutes = tabuSearch(data, initRoutes, maxIterations, maxTime, maxEvaluations); 
    auto mainEndTime = chrono::steady_clock::now();
    double executionTime = chrono::duration_cast<chrono::duration<double>>(mainEndTime - mainStartTime).count();

    if (bestRoutes.empty() && data.customers.size() > 1) {
        cout << "Final solution is empty." << endl;
    } else if (!bestRoutes.empty() && isSolutionFeasible(bestRoutes, data)) {
        cout << "Final solution is FEASIBLE." << endl;
    } else if (!bestRoutes.empty()) {
        cout << " Final solution is INFEASIBLE!" << endl;}
    int routeNumber = 1;
    for (size_t i = 0; i < bestRoutes.size(); ++i) {
        if (!bestRoutes[i].customerIds.empty()) {
            cout << "Route " << routeNumber++ << ":";
            for (int custId : bestRoutes[i].customerIds) {
                cout << " " << custId;}
            cout << endl;}}
    if (bestRoutes.empty() && data.customers.size() > 1) {
        cout << "(No routes in final solution)" << endl;}
    auto [finalVehicles_check, finalDistance_check] = objectiveFunction(bestRoutes);
    cout << "Vehicles: " << finalVehicles_check << endl;
    cout << "Distance: " << fixed << setprecision(2) << finalDistance_check << endl;
    cout << "Total Execution Time: " << fixed << setprecision(2) << executionTime << " seconds" << endl;
    return 0;
}