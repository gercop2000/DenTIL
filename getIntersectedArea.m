function A=getIntersectedArea(pointsShape1,pointsShape2)

polyarray1 = polyshape(pointsShape1);
polyarray2 = polyshape(pointsShape2);
polyout = intersect([polyarray1 polyarray2]);
A=polyarea(polyout.Vertices(:,1),polyout.Vertices(:,2));

end