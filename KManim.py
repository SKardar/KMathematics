# Manim:
from manim import * # type: ignore
# ===================================== Section 1 ======================================
# ****************************** Rate Function Generator *******************************
def simple_cubic_bezier(x0,y0,x1,y1):
    return bezier([
        np.array([0 , 0 , 0]),
        np.array([x0, y0, 0]),
        np.array([x1, y1, 0]),
        np.array([1 , 1 , 0]),
    ])   # type: ignore

def line_between_points(t, x0, y0, x1, y1):
    return (y1 - y0) * (t - x0) / (x1 - x0) + y0 # Linear Function : y-y0 = m(x-x0)

def closure_func(func, partitions=200):
    dt = 1 / partitions
    t_range = np.arange(0, 1+dt, dt)
    # x = f(x)[0] , y = f(x)[1]
    x_range = [func(t)[0] for t in t_range]
    y_range = [func(t)[1] for t in t_range]
    x_steps = [
        (x_range[i], x_range[i+1])
        for i in range(len(x_range)-1)
    ]
    y_steps = [
        (y_range[i], y_range[i+1])
        for i in range(len(y_range)-1)
    ]
    def f(t):
        if t <= 0:
            return 0
        elif t >=1:
            return 1
        for i in range(len(x_range)):
            y0 , y1 = y_steps[i]
            x0 , x1 = x_steps[i]
            if x0 <= t <= x1:
                return line_between_points(t, x0, y0, x1, y1)
    return f
# ++++++++++++++++++++++++ Section 1 Important Functions +++++++++++++++++++++++++++++++
def cubic_bezier(*args):
    return closure_func(simple_cubic_bezier(*args), partitions=200)

def smooth_cubic_bezier(*args):
    return lambda t: cubic_bezier(*args)(smooth(t))
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# **************************************************************************************
# ================================== End of Section 1 ==================================

# ===================================== Section 2 ======================================
# ***************** Make one animation start in the middle of another ******************
def custom_time(t,partitions,start,end,func):
    duration = end - start
    fragment_time = 1 / partitions
    start_time = start * fragment_time
    end_time = end * fragment_time
    duration_time = duration * fragment_time
    def fix_time(x):
        return (x - start_time) / duration_time
    if t < start_time: 
        return func(fix_time(start_time))
    elif start_time <= t < end_time:
        return func(fix_time(t))
    else:
        return func(fix_time(end_time))

# ++++++++++++++++++++++++ Section 2 Important Functions +++++++++++++++++++++++++++++++
def Custom(partitions,start,end,func=smooth):
    return lambda t: custom_time(t,partitions,start,end,func)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# **************************************************************************************
# ================================== End of Section 2 ==================================

# ===================================== Section 3 ======================================
# **************************** Defining Some Rate Functions ****************************
# Important Note:
# You can define additional parameteres as long as they have default values.
# Define 'amp' value here before using in your codes.
def parabola(t, amp = 1):
    return (1 - (2 * t - 1) ** 2) * amp

def multiline(t, amp = 1):
    if t < 0.5:
        return 2 * t * amp
    else:
        return (-2 * t + 2) * amp


# **************************************************************************************
# ================================== End of Section 3 ==================================