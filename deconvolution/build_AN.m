function AN = build_AN(Kq,q2)
    global MSpl L R N
    AN=-MSpl\(L+Kq);
    %AN = -MSpl\(L+R+Kq); % different boundary condition here leads to different AN.
    assert(all(size(AN)==[N+1,N+1]));
end