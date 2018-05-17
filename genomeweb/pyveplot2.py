"""
The creation of a Hive plot with Pyveplot is a straightforward process.
A Hive plot consists of:

    * radialy distributed linear axes
    * nodes along those axes
    * conections among those nodes

Pyveplots provides the corresponding objects: a *Hiveplot* class which
contains an arbitrary number of *Axis* objects which in turn contain an
arbitrary number of *Node* objects, and a method to connect them.

Change Notes:

    - Hiveplot, Axis and Node all inherit from Drawing class
      -> increase user power and decrease verbosity
      
         - User can parse Drawing options on init, e.g. size=('x%','y%')
           See svgwrite.drawing.Drawing() for more information
    
         - Allows addition of axis using Hiveplot.add(axis) rather than
           defining self.axes list
       
         - User no longer has to use Hiveplot.dwg.attribute to access
           Drawing elements
    
    - Removed Hiveplot.draw_axes() method; containg code was moved to
      Hiveplot.save()
      
    - Removed Axis.draw(); axis "draw" code moved to init and nodes are
      added to the axis when they are defined using Axis.add_node()
    
    - Removed Axis.getDwg(); same thing is achieved on calling Axis()
      once initiated
    
    - Aixs.angle is now property
      -> removes requirement for function call and can simply get 
         attribute using Axis.angle
    
    - Removed Node.getDwg() - just returns self
    
    - Removed example from Doc String. Example given did not work, see
      end for working example
      -> May be worth providing test suite and examples
    
    - Gave classes a simple __repr__
    
    - Applied PEP8 and some stylistic tweaks (subjective)
"""

from svgwrite import cm, mm, Drawing
from math import sin, cos, atan2, degrees, radians, tan, sqrt

class Hiveplot(Drawing):
    """
    Base class for a Hive plot.    
    """

    def __init__(self, *a, **kw):
        super(Hiveplot, self).__init__(*a, **kw)
        self.axes  = []
    
    def __repr__(self):
        return '<%s object(%d Axes, %d Nodes)>' % (
            self.__class__.__name__,
            len(self.axes),
            sum(len(ax.nodes) for ax in self.axes))
    
    def connect(
            self,
            axis0, n0_index, source_angle,
            axis1, n1_index, target_angle,
            curved=True,
            **kwargs):
        """
        Draw edges as Bezier curves.
        
        Start and end points map to the coordinates of the given nodes
        which in turn are set when adding nodes to an axis with the
        Axis.add_node() method, by using the placement information of
        the axis and a specified offset from its start point.
        
        Control points are set at the same distance from the start (or
        end) point of an axis as their corresponding nodes, but along
        an invisibleaxis that shares its origin but diverges by a given
        angle.
        
        Parameters
        ----------
        axis0        : source Axis object
        n0_index     : key of source node in nodes dictionary of axis0
        source_angle : angle of departure for invisible axis that 
                       diverges from axis0 and holds first control points
        axis1        : target Axis object
        n1_index     : key of target node in nodes dictionary of axis1
        target_angle : angle of departure for invisible axis that
                       diverges from axis1 and holds second control points
        kwargs       : extra SVG attributes for path element, optional
                       Set or change attributes using key=value
        """ 
        
        n0 = axis0.nodes[n0_index]
        n1 = axis1.nodes[n1_index]
        
        if curved:
            # source
            pth = self.path(d="M %d %d" % (n0.x, n0.y), fill='none', **kwargs)
            
            # compute source control point
            alfa = axis0.angle + source_angle
            length = sqrt(
                (n0.x - axis0.start[0])**2 + (n0.y - axis0.start[1])**2)
            x = axis0.start[0] + length * cos(alfa);
            y = axis0.start[1] + length * sin(alfa);
            
            # first control point in path
            pth.push("C %.f %.f" % (x, y))
            
            # compute target control point
            alfa = axis1.angle + target_angle
            length = sqrt(
                (n1.x - axis1.start[0])**2 + (n1.y - axis1.start[1])**2)
            x = axis1.start[0] + length * cos(alfa);
            y = axis1.start[1] + length * sin(alfa);
            
            # second control point in path
            pth.push("%.f %.f" % (x, y))
            # target
            pth.push("%.f %.f" % (n1.x, n1.y))
            self.add(pth)
        else:
            self.add(self.line(
                (n0.x, n0.y),
                (n1.x, n1.y),
                **kwargs))
      
    def save(self):
        # add any remaining axis elements before saving
        for axis in self.axes:
            self.add(axis())
        super(Hiveplot, self).save()
        
class Axis(Drawing):
    """
    Initialize Axis object with start, end positions and optional SVG
    attributes
    
    Parameters
    ----------
    start  : (x, y) start point of the axis
    end    : (x1, y1) end point of the axis
    kwargs : extra SVG attributes for line element, optional
             Set or change attributes using key=value
    """
    
    def __init__(self, start=(0, 0), end=(0, 0), **kwargs):
        super(Axis, self).__init__()
        self.start = start
        self.end = end
        self.nodes = []
        self.kw = kwargs
    
    def __call__(self):
        # finalize object
        self.add(
            self.line(
                start=self.start,
                end=self.end,
                **self.kw))
        return self
    
    def __repr__(self):
        return '<%s object(%d Nodes)>' % (
            self.__class__.__name__,
            len(self.nodes))
    
    def __len__(self):
        return sqrt(
            (self.end[0] - self.start[0])**2 +
            (self.end[1] - self.start[1])**2)
    
    def add_node(self, node, offset):
        """
        Add a Node object to nodes dictionary, calculating its
        coordinates using offset
        
        Parameters
        ----------
        node   : a Node object
        offset : float 
                 number between 0 and 1 that sets the distance
                 from the start point at which the node will be placed
        """
        width = self.end[0] - self.start[0]
        height = self.end[1] - self.start[1]        
        node.x = self.start[0] + (width * offset)
        node.y = self.start[1] + (height * offset)
        
        if node.draw_label:
            node.add_label(self.start[0])
        
        self.nodes.append(node)
        self.add(node)
    
    @property
    def angle(self):
        xDiff = self.end[0] - self.start[0]
        yDiff = self.end[1] - self.start[1]
        return atan2(yDiff, xDiff)

class Node(Drawing):
    """
    Base class for Node objects.
    
    Holds coordinates for node placement and a svgwrite.Drawing()
    object.
    
    Parameters
    ----------
    ID: a unique key for the nodes dict of an axis.  
    """
    def __init__(self, ID='', draw_label=False, *a, **kw):
        super(Node, self).__init__(*a, **kw)
        self.ID = ID
        self.x = 0
        self.y = 0
        self.r = 1.5
        self.draw_label = draw_label
        
    def add_label(self, ox=0):
        self.add(self.text(
            self.ID,
            insert=(self.x, self.y),
            style='font-size:11px; font-family:Arial',
            text_anchor='end' if self.x < ox else 'start'))

if __name__ == '__main__':
    # A short example
    import networkx
    import random
    # a network
    g = networkx.barabasi_albert_graph(50, 2)
    
    # hiveplot object
    h = Hiveplot('short_example.svg', size=('200%, 200%'))
               # start     end
    h.axes = [
        Axis((200,200), (200,100), stroke="grey"),
        Axis((200,200), (300,300), stroke="blue"),
        Axis((200,200), (10,310), stroke="black")]
    
    # randomly distribute nodes in axes
    for n in g.nodes():
        random.choice(h.axes).add_node(Node(n), random.random())
    
    for e in g.edges():
        if (e[0] in h.axes[0].nodes) and (e[1] in h.axes[1].nodes):
            # edges from axis0 to axis1    
            h.connect(
                h.axes[0], e[0], 45,
                h.axes[1], e[1], -45,
                stroke_width='0.34', stroke_opacity='0.4', stroke='purple')
        elif (e[0] in h.axes[0].nodes) and (e[1] in h.axes[2].nodes):
            # edges from axis0 to axis2
            h.connect(
                h.axes[0], e[0], -45,
                h.axes[2], e[1], 45,
                stroke_width='0.34', stroke_opacity='0.4', stroke='red')
        elif (e[0] in h.axes[1].nodes) and (e[1] in h.axes[2].nodes):
            # edges from axis1 to axis2
            h.connect(
                h.axes[1], e[0], 15,
                h.axes[2], e[1], -15,
                stroke_width='0.34', stroke_opacity='0.4', stroke='magenta')
    h.save()
        
