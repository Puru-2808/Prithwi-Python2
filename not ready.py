import datetime
a=datetime.datetime.now()

print("Hi, I'm Purusottam's AI...I'm under construction")

n=str(input("What's your name: "))

if a.hour<12:
    print("Good morning",n,",how can i help you?")
elif a.hour==12:
    print("Good noon",n,",how can i help you?")
elif a.hour>12 and a.hour<18:
    print("Good afternoon",n,",how can i help you?")
elif a.hour>=18:
    print("Good evening",n,",how can i help you?")

b=str(input("lets do some simple arithmatic math, how is it...ready?(yes/no): "))
if b=="yes":
    while True:
         op=str(input("Enter operation type(+,-,×,÷): "))
         if op== '+':
              n1=float(input("Enter 1st number: "))
              n2=float(input("Enter 2nd number: "))
              s=n1+n2
              print("Result=",s)
         if op== '-':
              n1=float(input("Enter 1st number: "))
              n2=float(input("Enter 2nd number: "))
              s=n1-n2
              print("Result=",s)
         if op=='×':
             n1=float(input("Enter 1st number: "))
             n2=float(input("Enter 2nd number: "))
             s=n1*n2
             print("Result=",s)
         if op=='÷':
             n1=float(input("Enter 1st number: "))
             n2=float(input("Enter 2nd number: "))
             s=n1/n2
             print("Result=",s)
         a=str(input("Do you want to continue(yes/no)?: "))
         if a== "no":
            break
    print(".............Calculator is closing...............")
else:
    print("Okay, Actually i was in mathematics mood...but if you don't want to do this then i am closing...bye")   