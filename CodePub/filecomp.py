import ast

with open("exmat/a.23.txt", "r") as napov_ja:
    with open("mat/a.txt", "r") as ja:
        lines_napov = napov_ja.readlines()
        lines_ja = ja.readlines()
        
        for i in range(len(lines_napov)):
            if ast.literal_eval(lines_napov[i]) != ast.literal_eval(lines_ja[i]):
                print(i + 1)
