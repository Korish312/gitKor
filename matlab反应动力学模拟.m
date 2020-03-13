a=input('what is k1?\n');
b=input('what is k-1?\n');
c=input('what is k2?\n');
d=input('what is k-2?\n');

tspan=0:1000:5000;
y0=[1;1;0;0];
[t,y]=ode45(@(t,y) odefcn(t,y,a,b,c,d), tspan, y0);
plot(t,y(:,1),'r-',t,y(:,2),'g-',t,y(:,3),'b-',t,y(:,4),'k-')

function dy=odefcn(~,y,a,b,c,d)
dy=zeros(4,1);
dy(1)=-a*y(1)*y(2)+b*y(3);
dy(2)=-a*y(1)*y(2)+b*y(3);
dy(3)=a*y(1)*y(2)-c*y(3)-b*y(3)+d*y(4);
dy(4)=c*y(3)-d*y(4);
end
