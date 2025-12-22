def f(a,b,d,password):
    for i in range(len(a)):
         if a[i]==d and b[i]==password:
             return i
    return 0

a=["a","b"]
b=[1,2]
d=str(input("enter name: "))
password=int(input("enter password: "))

p=f(a,b,d,password)

if p!=0:
    print("success")
else:
    print("error")