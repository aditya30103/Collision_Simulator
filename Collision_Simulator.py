class PriorityQueueBase:
# Abstract base class for a priority queue
    # Class defined to store values together
    class _Node:
        __slots__ = '_time' , '_object' , '_position', '_index'
        def __init__ (self, t, i, x, ind):
            self._time = t
            self._object = i
            self._position = x
            self._index = ind

        # OPERATOR OVERLOADING
        def __lt__ (self, other): # compare items based on their keys
            if (self._time == 'No_Collision'): return False
            if (other._time == 'No_Collision'): return True
            if (self._time < other._time): return True
            if (self._time == other._time and self._object < other._object): return True
            if (self._time == other._time and self._object == other._object and self._index < other._index): return True
            return False

class HeapPriorityQueue(PriorityQueueBase): # Base Class Defines Node
# Using Min-Heap to Implement Priority Queue
#------------------------------ Non-Public Functions ------------------------------
    def _parent(self, j):
        return (j-1) // 2

    def _left(self, j):
        return 2*j+1

    def _right(self, j):
        return 2*j+2

    def _has_left(self, j):
        return self._left(j) < len(self._data) # index beyond end of list?

    def _has_right(self, j):
        return self._right(j) < len(self._data) # index beyond end of list?

    def _swap(self, i, j):
        self._data[i], self._data[j] = self._data[j], self._data[i]
        self._data[i]._index = i
        self._data[j]._index = j
        self._refer[self._data[i]._object]=i
        self._refer[self._data[j]._object]=j

    def _upheap(self, j):
        parent = self._parent(j)
        if (j > 0 and self._data[j] < self._data[parent]):
            self._swap(j, parent)
            self._upheap(parent) # recur at position of parent

    def _downheap(self, j):
        if (self._has_left(j)):
            left = self._left(j)
            small_child = left # although right may be smaller
            if self._has_right(j):
                right = self._right(j)
                if self._data[right] < self._data[left]:
                    small_child = right
            if self._data[small_child] < self._data[j]:
                self._swap(j, small_child)
                self._downheap(small_child) # recur at position of small child
    
    def _bubble(self, j):
        if j > 0 and self._data[j] < self._data[self._parent(j)]:
            self._upheap(j)
        else:
            self._downheap(j)

#------------------------------ Public Functions ------------------------------
    def __init__ (self, contents=()): # Constructor
        self._data = [] 
        self._refer = []
        for val in contents:
            obj = self._Node(val[0], val[1], val[2], len(self._data))
            self._refer.append(len(self._data))
            self._data.append(obj)            
        if len(self._data) > 1:
            self._heapify()

    def _heapify(self):
        start = self._parent(len(self) - 1) # start at PARENT of last leaf
        for j in range(start, -1, -1): # going to and including the root
            self._downheap(j)

    def __len__ (self):
        return len(self._data)

    def min(self):
        item = self._data[0]
        return (item._time, item._object, item._position)

    def update_heap(self, new_tix):
        for i in range(0,len(new_tix)):
            obj_no = new_tix[i][1]
            obj_ref = self._refer[obj_no]
            new_node = self._Node(new_tix[i][0], new_tix[i][1], new_tix[i][2], obj_ref)
            self._data[obj_ref] = new_node
            self._bubble(obj_ref)

#---------------------------------------- End Of Data Structure ------------------------------------
#---------------------------------------- Helper Functions For ListCollisions ----------------------

def calculator(M, v, obj_ind):
    m1 = M[obj_ind]
    m2 = M[obj_ind+1]

    u1 = v[obj_ind]
    u2 = v[obj_ind+1]

    Num_v1 = (m1-m2)*u1 + (2*m2*u2)
    Num_v2 = (2*m1*u1) - (m1-m2)*u2
    Den = (m1+m2)

    v1 = Num_v1/Den
    v2 = Num_v2/Den
    return [v1, v2]

def update_x_lastcolt(x, v, lastcolt, time_counter, pos_col, obj_ind):
    if (obj_ind-1 > -1): 
        x[obj_ind-1] = x[obj_ind-1] + v[obj_ind-1]*(time_counter-lastcolt[obj_ind-1])
        lastcolt[obj_ind-1] = time_counter
    x[obj_ind] = pos_col
    lastcolt[obj_ind] = time_counter
    x[obj_ind+1] = pos_col
    lastcolt[obj_ind+1] = time_counter
    if (obj_ind+2 <= len(x)-1): 
        x[obj_ind+2] = x[obj_ind+2] + v[obj_ind+2]*(time_counter-lastcolt[obj_ind+2])
        lastcolt[obj_ind+2] = time_counter

def roundoff(tup):
    a = [0,0,0]
    for i in range(0,3):
        a[i] = round(tup[i], 4)
    rtup = tuple(a)
    return rtup

def new_tix_generator(M, x, v, time_counter, obj_ind):
    new_tix = []
    for i in range(-1,2):
        if (obj_ind+i > -1 and obj_ind+i+1 <= len(M)-1):
            dis = x[obj_ind+i+1] - x[obj_ind+i]
            rvel = v[obj_ind+i] - v[obj_ind+i+1]
            if (rvel == 0.00): time = 'No_Collision'
            else:
                time = dis/rvel + time_counter
                if (time <= time_counter): time = 'No_Collision'                
            new_tix.append((time, obj_ind+i, x[obj_ind+i]))
    return new_tix

def listCollisions(M, x, v, m, T): #------ Time Complexity: O(n + mlogn) --------------------
    
    tix = [] 
    lastcolt = [0]*len(M)
    for i in range(0,len(M)-1): #-------- Constructing tix tuples in O(n) -------------------
        dis = x[i+1] - x[i]
        rvel = v[i] - v[i+1]
        if (rvel == 0.00): time = 'No_Collision'
        else:
            time = (dis)/(rvel)
            if (time <= 0.00): time = 'No_Collision'
        tix.append((time, i, x[i]))

    A = HeapPriorityQueue(tix) #--------- Heapifying tix in O(n) ----------------------------

    Collisions = []
    time_counter = 0.00
    
    for i in range(0,m): #--------------- Iterating Over Number Of Max_Collisions -----------
        
        col = A.min()
        if (col[0] == 'No_Collision'): return Collisions  #---------------- No more Collisions Left -------
        if (col[0] > time_counter): time_counter = col[0]
        if (time_counter > T) : return Collisions #------------------------ Time Limit Reached ------------
        
        obj_ind = col[1]
        f_col = (col[0], obj_ind, x[obj_ind] + v[obj_ind]*(time_counter-lastcolt[obj_ind]))    
        pos_col = f_col[2]

        rf_col = roundoff(f_col) # rounding off value to be printed to 4 decimal places
        Collisions.append(rf_col) # appending in list Collisions
     
        calc = calculator(M, v, obj_ind) # calculating new velocities of object which collided

        v[obj_ind] = calc[0]
        v[obj_ind+1] = calc[1]
        
        update_x_lastcolt(x, v, lastcolt, time_counter, pos_col, obj_ind) # updating list x and lastcolt
        new_tix = new_tix_generator(M, x, v, time_counter, obj_ind) # generating new tix tuples for futher collisions   
        A.update_heap(new_tix) # heap updated for fresh iteration
    return Collisions