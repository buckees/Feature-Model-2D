"""
Contain Misc. Functions
"""

# A utility function to calculate area  
# of triangle formed by (x1, y1),  
# (x2, y2) and (x3, y3) 
  
def trgl_area(x1, y1, x2, y2, x3, y3): 
  
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1)  
                + x3 * (y1 - y2)) / 2.0) 
  
  
# A function to check whether point P(x, y) 
# lies inside the triangle formed by  
# A(x1, y1), B(x2, y2) and C(x3, y3)  
def isInsideTrgl(x1, y1, x2, y2, x3, y3, x, y): 
  
    # Calculate area of triangle ABC 
    A = trgl_area(x1, y1, x2, y2, x3, y3) 
  
    # Calculate area of triangle PBC  
    A1 = trgl_area(x, y, x2, y2, x3, y3) 
      
    # Calculate area of triangle PAC  
    A2 = trgl_area(x1, y1, x, y, x3, y3) 
      
    # Calculate area of triangle PAB  
    A3 = trgl_area(x1, y1, x2, y2, x, y) 
      
    # Check if sum of A1, A2 and A3  
    # is same as A 
    if(A == A1 + A2 + A3): 
        return True
    else: 
        return False

if __name__ == '__main__':  
    # Driver program to test above function 
    # Let us check whether the point P(10, 15) 
    # lies inside the triangle formed by  
    tempx, tempz = 25.0, 400.0
    
    if (isInsideTrgl(30.0, 350.0, 30.0, 450.0, 0.0, 450.0, tempx, tempz)): 
        print('Inside') 
    else: 
        print('Not Inside')  
      
    # This code is contributed by Danish Raza 
    
