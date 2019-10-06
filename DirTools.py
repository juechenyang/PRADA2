"""
created by Juechen Yang at 10/4/18

"""
import os

def check_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def make_dir(path):
    os.makedirs(path)

