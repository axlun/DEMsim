import sys
import os



if __name__ == '__main__':
    argument_string = sys.argv[1]

    print('Hello world')

    print(argument_string)
    print(os.listdir(argument_string))