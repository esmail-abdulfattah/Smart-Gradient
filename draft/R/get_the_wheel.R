

#0.5: from 0 to pi/2 and 2 is for full wheel
#It: getting the value of x at Iteration "It" using Global
It = 4
resrot <- rotate_the_Wheel(It, Global, wheel_parts = 2)
plot_rotated_track(resrot$errors,"",Global$new.trace[It],Global$default.trace[It])
get_the_wheel(resrot$vects,resrot$errors)

