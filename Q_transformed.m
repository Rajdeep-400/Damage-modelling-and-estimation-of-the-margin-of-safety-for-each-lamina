% Function to calculate reduced stiffness matrix (plane stress) for a given theta
function Qt = Q_transformed(Q, m, n)
    Tg = [m.^2, n.^2, 2.*m.*n; n.^2, m.^2, -2.*m.*n; -m.*n, m.*n, (m.^2-n.^2)]; %transformation matrix for stress matrix
    Te = [m.^2, n.^2, m.*n; n.^2, m.^2, -m.*n; -2.*m.*n, 2.*m.*n, (m.^2-n.^2)]; %transformation matrix for strain matrix
    Qt = inv(Tg)*Q*Te;
end