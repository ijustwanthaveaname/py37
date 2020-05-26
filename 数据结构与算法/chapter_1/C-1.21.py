textbox=[]
try:
    while True:
        text=input('please input message: ')
        textbox.append(text)
except EOFError:
    textbox.reverse()
    for i in textbox:
        print(i)


