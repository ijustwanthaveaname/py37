import os
def dist_usage(path):
    total=os.path.getsize(path)
    if os.path.isdir(path):
        for filename in os.listdir(path):
            childpath=os.path.join(path,filename)
            total+=dist_usage(childpath)
    print('{0:<7}'.format(total),path)
    return total