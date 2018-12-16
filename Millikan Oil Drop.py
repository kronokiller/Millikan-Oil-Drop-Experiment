import os
import math
import statistics
import matplotlib.pyplot as plt
from scipy.stats import linregress
from numpy import linspace

# Open the file containing the data assuming the file is in the same directory
os.chdir(os.path.dirname(__file__))
file = open('Millikan Oil Drop.csv')

# Create a list to contain each data point which consists of whether the drop is fallinig or rising,
# the resistance of the thermistor, the time that it took, and the number of major reticle markings
# that were traversed
data = []
for line in file:
    # Get a list of the 4 items separated by commas
    list = line.split(',')
    # Remove the new line character from the last value
    list[-1] = list[-1][:-1]
    # 1s and 0s were used to record whether the drop was rising or falling
    if list[0] == '1':
        list[0] = 'fall'
    if list[0] == '0':
        list[0] = 'rise'
    # The thermistor resistance was a float
    list[1] = float(list[1])
    # The time was recorded as seconds.frames since there were 30 frames per second, the frames were
    # divided by frame rate and added to the seconds.
    list[2] = int(float(list[2])) + float(list[2]) % 1 / 30
    # The amount of major reticle markings that were passed was given as an integer between 1 and 4
    list[3] = int(list[3])
    # Each list will be converted into a tuple to represent each data point
    data.append(tuple(list))

# The conversions values for the thermistor are given in ohms for temperatures in 1 K intervals starting at 283.15 K
thermistorTable = [3.239, 3.118, 3.004, 2.897, 2.795, 2.7, 2.61, 2.526, 2.446, 2.371, 2.3, 2.233, 2.169, 2.11, 2.053, 2, 1.95, 1.902, 1.857, 1.815, 1.774, 1.736, 1.7, 1.666, 1.634, 1.603, 1.574, 1.547, 1.521, 1.496]


# Calculates the temperature based on the resistance of the thermistor
def T(omega):
    i = 0
    # Iterates over the pairs of resistances in the conversion table
    while not(omega <= thermistorTable[i] and omega > thermistorTable[i + 1]):
        # If the resistance isn't between the pair, check the next pair
        i += 1
    # The index of the lower temperature in the pair is added to the temperature of the first index
    # The temperature is approximated as a linear function and the fraction of the way to the next index is added to the temperature
    return i + 283.15 + (omega - thermistorTable[i]) / (thermistorTable[i + 1] - thermistorTable[i])

# Calculates the viscosity of the air
def nu(T):
    # The viscosity of the air is a linear function of temperature between the points (288.15, 1.8 * 10^(-5)) and (305.15, , 1.881 * 10^(-5))
    return (T - 288.15) / 17 * 0.081 * 10 ** (-5) + 1.8 * 10 ** (-5)

# Calculates the pressure
def p(T):
    # Commented equation may not be valid because pressure in the building is likely different from outside pressure and
    # various other factors such as not knowing the environmental conditions
    #return Po * math.e ** (-mu * g * h / R / T)
    return Po

# Calculates the average fall velocity
def vf(data):
    # Creates a list of fall velocity calculations
    vList = []
    # Iterates over all of the data points in the data
    for point in data:
        # only calculates fall velocity for falling drops
        if point[0] == 'fall':
            # calculates the speed as the distance per reticle marking times the amount of reticle markings divided by the tiem
            vList.append(l * point[3] / point[2])
    # Calculates the average of the list
    return sum(vList) / len(vList)

# Calculates the standard deviation of the average fall velocity
def vfStDev(data):
    # Creates a list of fall velocity calculations
    vList = []
    # Iterates over all of the data points in the data
    for point in data:
        # only calculates fall velocity for falling drops
        if point[0] == 'fall':
            # calculates the speed as the distance per reticle marking times the amount of reticle markings divided by the tiem
            vList.append(l * point[3] / point[2])
    # Calculates the standard deviation of the list
    return statistics.stdev(vList)

# Calculates the rise velocity using the time and number of reticle markings passed
def vr(t, n):
    return l * n / t

# Calulates the excess charge of the drop using the resistance, the time, the number of reticals markiings passed,
# other global quantities, and the data used for the fall velocity
def Q(omega, t, n, data):
    return 4 * math.pi / 3 * rho * g * d * (vf(data) + vr(t, n)) / V / vf(data) * (((b / 2 / p(T(omega))) ** 2 + 9 * nu(T(omega)) * vf(data) / 2 / rho / g) ** (1 / 2) - b / 2 / p(T(omega))) ** 3

# The density of the oil
rho = 886

# The distance between major reticle markings
l = .0005

# The approximate sea level pressure
Po = 101325

# The universal gas constant for air
R = 8.31432

# The altitude of Charleston, Il
h = 211

# The molar mass of Earth's air
mu = 0.0289644

# The distance between the plates
d = .00713

# A constant derived by Millikan
b = .0082

# The acceleration due to gravity
g = 9.8

# The voltage used between the plates
V = 500

print(vf(data), vfStDev(data))

# Get the excess charge and velocity of the oil drop for every rising drop
QList = []
for point in data:
    if point[0] == 'rise':
        QList.append([Q(point[1], point[2], point[3], data), vr(point[2], point[3])])

# For each item in the list of excess charges, approximate the number of excess electron using the given value
# and calculate the charge of a single electron by dividing the total excess charge by the number of excess electrons
for i in range(len(QList)):
    # The approximate number of excess electrons
    N = QList[i][0] / (1.6 * 10 ** (-19))
    # The charge of the electron
    q = QList[i][0] / N
    QList[i] = [QList[i][0], QList[i][1], N, q]

# Plot each charge as a point along the line y = x
plot = plt.scatter([item[1] for item in QList], [item[1] for item in QList])
plt.axis(ymax = 0.001, ymin = 0, xmin = 0, xmax = 0.001)
plt.show()

# Group each according to visuals seen in plot above
groupings = [[] for i in range(12)]
for i in range(len(QList)):
    if QList[i][1] < .000034 and QList[i][1] > 0.000026:
        groupings[0].append(QList[i])
    elif QList[i][1] < .0000940 and QList[i][1] > 0.0000905:
        groupings[1].append(QList[i])
    elif QList[i][1] < .00016 and QList[i][1] > 0.00014:
        groupings[2].append(QList[i])
    elif QList[i][1] < .00026 and QList[i][1] > 0.0002:
        groupings[3].append(QList[i])
    elif QList[i][1] < .00031 and QList[i][1] > 0.00028:
        groupings[4].append(QList[i])
    elif QList[i][1] < .00034 and QList[i][1] > 0.00032:
        groupings[5].append(QList[i])
    elif QList[i][1] < .00038 and QList[i][1] > 0.00037:
        groupings[6].append(QList[i])
    elif QList[i][1] < .000401 and QList[i][1] > 0.000398:
        groupings[7].append(QList[i])
    elif QList[i][1] < .000502 and QList[i][1] > 0.000496:
        groupings[8].append(QList[i])
    elif QList[i][1] < .000668 and QList[i][1] > 0.000664:
        groupings[9].append(QList[i])
    elif QList[i][1] < .000750 and QList[i][1] > 0.000748:
        groupings[10].append(QList[i])
    elif QList[i][1] < .001 and QList[i][1] > 0.00099:
        groupings[11].append(QList[i])

# Locate any groups that have less than 5 points
popList = []
for l in groupings:
    if len(l) < 5:
        popList.append(groupings.index(l))

# Remove these groups of data points
popList.reverse()
for i in popList:
    groupings.pop(i)

# Calculate the mean and standard deviation of each group
calculations = []
for l in groupings:
    points = len(l)
    mean = sum([i[0] for i in l]) / len(l)
    standardDev = statistics.stdev([i[0] for i in l])
    calculations.append([mean, standardDev, points])

# Sort the calculated charges in ascending order
calculations = sorted(calculations)

# Display the calculated charges, difference between each, and uncertainties
print('')
print(calculations[0][0], '+/-', calculations[0][1])
for i in range(len(calculations) - 1):
    difference = calculations[i + 1][0] - calculations[i][0]
    differenceStDev = (calculations[i + 1][1] ** 2 + calculations[i][1] ** 2) ** (1 / 2)
    print('\n       ', difference, '+/-', differenceStDev)
    print('')
    print(calculations[i + 1][0], '+/-', calculations[i + 1][1])

# Looking at data the amounts of excess electrons appear to be these
numbers = [1, 4, 5, 7, 8]

# Plot the charges vs the number of excess electros and plot a linear regression
# linregress calculates the slope, y-intercept, r-value, p-value, and uncertainty in the slope
linearRegression = linregress(numbers, [i[0] for i in calculations])
plt.scatter(numbers, [i[0] for i in calculations])
plt.plot(linspace(0, 9, num = 200), linearRegression[1] + linearRegression[0] * linspace(0, 9, num = 200), 'r', label = 'q = ' + str(round(linearRegression[0] * 10 ** 22) / 10 ** 22) + ' * N - ' + str(round(-linearRegression[1] * 10 ** 22) / 10 ** 22))
plt.axis(ymax = 1.5 * 10 ** -18, ymin = 0, xmin = 0, xmax = 9)
plt.legend()
plt.show()

# print the slope with its error
print('\n')
print('q = m =', '(' + str(round(linearRegression[0] * 10 ** 21) / 100), '+/-', str(round(linearRegression[4] * 10 ** 21) / 100) + ')', '* 10E-19 C')