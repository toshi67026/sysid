function dx = f_c(x, u)
    global l1 l2 a1 a2 m1 m2 d1 d2 g phi1 phi2 phi3
    theta1 = x(1);
    theta2 = x(2);
    dtheta1 = x(3);
    dtheta2 = x(4);
    m2 = x(5);

    D = [phi1 + phi2 + 2 * phi3 * cos(theta2), phi2 + phi3 * cos(theta2);
        phi2 + phi3 * cos(theta2), phi2];

    h = (-phi3 * sin(theta2) * [dtheta2, dtheta1 + dtheta2; -dtheta1, 0] + [d1, 0; 0, d2]) * [dtheta1; dtheta2];

    g_q = g * [m1 * l1 * cos(theta1) + m2 * (a1 * cos(theta1) + l2 * cos(theta1 + theta2));
            m2 * l2 * cos(theta1 + theta2)];

    dx = [dtheta1; dtheta2;
        inv(D) * (u - h - g_q); 0];
end
