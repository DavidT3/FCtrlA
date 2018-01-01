import math
import sys


def r500(t, z):
    return 1.104 * math.pow(t/5.0, 0.57) * math.pow(scale(z), -0.5)


def scale(z):
    return 0.27*math.pow(1+z, 3) + 0.73


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        print('enter temp and redshift')
        sys.exit(0)

    temp = float(sys.argv[1])
    red = float(sys.argv[2])
    print(r500(temp, red) * 1000)



