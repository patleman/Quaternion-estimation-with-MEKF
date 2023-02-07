function result = Attiesticontrol( tstart, tend, dt)
    % Physical constants.
    g = 9.81;
    I = diag([5e-3, 5e-3, 10e-3]);
    % Simulation times, in seconds.
    if nargin < 4
        tstart = 0;
        tend = 12;
        dt = 0.001;
    end
    ts = tstart:dt:tend;
    % Number of points in the simulation.
    N = numel(ts);
    % Output values, recorded as the simulation runs.
    %thetaout = zeros(3, N);
    omegaout = zeros(3, N);
    inputMout = zeros(3, N);
    dentaout=zeros(3,N);
    quatout=zeros(4,N);
    errorout=zeros(4,N);

% Initial system state.
      Pentaq=[0.1 0 0;0 0.1 0;0 0 0.1];
      Rw = eye(3,3)*1e-6;
      Rv = eye(3,3)*1e-6;
      quat = [0; 0; 0; 1];
      omega=[0;0;0];
      P = 1e1*eye(6,6);
      Qw = eye(3,3)*1e-1;
      Qv = eye(3,3)*1e-2; 
      denta=[10^-7;10^-7;10^-7];
      venta=[10^-3;10^-3;10^-3];
      ind = 0;
      M=[0;0;0];
      errsum=0;
      errsum0=0;
      errsum1=0;  
    for t = ts
        %%%%%%%
        %%% to simulated values derived through equations, we will add
        %%% noise to them and then assume it to be measurement values
        %%%these measurement values go through MEKF and will finally give
        %%% quaternion and omega estimate that will be further taken into
        %%% consideration while computation.
        %%%  making measurement values ,we will be measuring omega and
        %%% a vector like gravity or magnetic field 
        Rot1 = QuatToRot(quat);
       [phies, thetes, psies] = RotToRPY_ZXY(Rot1); 
        Dz= [phies, thetes, psies]';
        qw=quat(1);
        q1=quat(2);
        q2=quat(3);
        q3=quat(4);
        quat=[q1;q2;q3;qw];
        truequat=quat;
        v=[0;0;-g];
         RT = [-quat(2)*quat(2)-quat(3)*quat(3)  quat(1)*quat(2)+quat(3)*quat(4)                                            quat(1)*quat(3)-quat(2)*quat(4);
             quat(1)*quat(2)-quat(3)*quat(4)  -quat(1)*quat(1)-quat(3)*quat(3)   quat(2)*quat(3)+quat(1)*quat(4);
             quat(1)*quat(3)+quat(2)*quat(4)   quat(2)*quat(3)-quat(1)*quat(4)  -quat(1)*quat(1)-quat(2)*quat(2)];
         RT = RT + RT + eye(3);   
        % we generate an additional vector acting on the system  
        nv=10^-7;
        rv=10^-7;
        rw=10^-7;
        qv = normrnd( 0 , nv , [3 1] );
        % we obtain the simulated measurements
        vm = RT*( qv + v ) + normrnd( 0 , rv , [3 1] );
        wm = omega + normrnd( 0 , rw , [3 1] );
        %%%%%% applying MEKF to above noisy measurements
        wnorm = norm(omega);
        if wnorm > 0
          wdt05 = 0.5*wnorm*dt;
           swdt = sin(wdt05)/wnorm;
           qw = [ omega*swdt; cos(wdt05) ];
        else
          qw = [0; 0; 0; 1];
        end
        qp = [quat(4)*qw(1:3) + qw(4)*quat(1:3) + cross( quat(1:3) , qw(1:3) );
              quat(4)*qw(4) - quat(1:3)'*qw(1:3) ];

% we compute the covariance matrix for the state prediction
         F = [-qw(2)*qw(2)-qw(3)*qw(3)   qw(1)*qw(2)+qw(3)*qw(4)   qw(1)*qw(3)-qw(2)*qw(4);
            qw(1)*qw(2)-qw(3)*qw(4)  -qw(1)*qw(1)-qw(3)*qw(3)   qw(2)*qw(3)+qw(1)*qw(4);
            qw(1)*qw(3)+qw(2)*qw(4)   qw(2)*qw(3)-qw(1)*qw(4)  -qw(1)*qw(1)-qw(2)*qw(2)];

         F = F + F + eye(3);
         Mj = [F           eye(3,3)*dt;
           zeros(3,3)  eye(3,3)   ];
      
         P = Mj*( P + [ zeros(3,3)   zeros(3,3);
                            zeros(3,3)   Qw*dt ] )*Mj';
      % we compute the measurement prediction
         RT = [ -qp(2)*qp(2)-qp(3)*qp(3)   qp(1)*qp(2)+qp(3)*qp(4)   qp(1)*qp(3)-qp(2)*qp(4);
                qp(1)*qp(2)-qp(3)*qp(4)  -qp(1)*qp(1)-qp(3)*qp(3)   qp(2)*qp(3)+qp(1)*qp(4);
                qp(1)*qp(3)+qp(2)*qp(4)   qp(2)*qp(3)-qp(1)*qp(4)  -qp(1)*qp(1)-qp(2)*qp(2) ];
         RT = RT + RT + eye(3,3);
      
         vp = RT*v;
      % we compute the covariance matrix for the measurement prediction
         F = [  0     -vp(3)   vp(2);
            vp(3)    0     -vp(1);
           -vp(2)   vp(1)    0   ];
      
         Mj = [F           zeros(3,3);
           zeros(3,3)  eye(3,3)  ];
      
         S = Mj*P*Mj' + [Qv + Rv  zeros(3,3);
                     zeros(3,3)       Rw  ];      
% now we can compute the Kalman gain
         K = P*Mj'/S;
      % and update the state in the chart
         dx = K*[vm-vp; wm-omega];
      % the updated point in the chart is mapped to a quaternion
         delta =fC2M( dx(1:3) );
         quat = [qp(4)*delta(1:3) + delta(4)*qp(1:3) + cross( qp(1:3) , delta(1:3) );
               qp(4)*delta(4) - qp(1:3)'*delta(1:3) ];
      % and the angular velocity is updated in the usual way
         omega = omega + dx(4:6);
      % the covariance matrix is updated in the chart centered in qp
          Mj = eye(6,6) - K*Mj;
          P = Mj*P;
       % finally we update the covariance matrix from the chart centered in
      %  qp quaternion, to the chart centered in the updated q quaternion
         F = chartUpdateMatrix( delta );
          Mj = [F           zeros(3,3);
           zeros(3,3)  eye(3,3)];
         P = Mj*P*Mj';
      % we avoid numerical instabilities
        quat = quat/norm(quat);
        estiquat=quat;
        P = 0.5*(P + P');
        error=truequat-estiquat;
        %%%%%%
        ind = ind + 1;
        qd(1)=0.5*sin(t);
        qd(2)= 0.5*cos(t);
        qdp(1)= 0.5*cos(t);
        qdp(2) = -0.5*sin(t);
        qw=quat(4);
        q1=quat(1); 
        q2=quat(2);
        q3=quat(3);
        quat=[qw;q1;q2;q3];
        Rot = QuatToRot(quat);
        [phi, thet, psi] = RotToRPY_ZXY(Rot);  %%  ZXY SYSTEM

        pol=[cos(thet) 0 -cos(phi)*sin(thet);
            0        1  sin(phi);
            sin(thet) 0  cos(phi)*cos(thet)];
        J_q=pol'*I*pol;
        mu=0;covwtilda=(10)^-7;wtilda=mvnrnd(mu,covwtilda);
        %%%%%%
        disturbance=[0.06*sin(t);0.06*sin(t);0.06*sin(t)]+[wtilda wtilda                     wtilda]';
        mu=0;covw=(10)^-10;w=mvnrnd(mu,covw);
        phi_des=qd(1);
        theta_des=qd(2);
        yawdot=0;
        yaw=1.5;
        p_des=qdp(1);q_des=qdp(2);
        p=omega(1);q=omega(2);r=omega(3);
        %%%%%%% disturbance oobserver
        L=5.0 ; %
        dventa=-L*(denta)+L*(cross(omega, I*omega)-M);
        venta=venta+dventa*dt;
        denta=venta+L*I*omega;   
     %%%%%% disturbance observer end PID STARTS
        errsum=errsum+yaw-psi;
        errsum1=errsum1+theta_des-thet;
        errsum0=errsum0+phi_des-phi;
        M = [120*(phi_des-phi)+1.289*(p_des-p)+ 0.002*errsum0;%%%%%%%%%moment
             120*(theta_des-thet)+1.289*(q_des-q)+ 0.002*errsum1;
             120*(yaw-psi)+1.289*(yawdot-r)+(.002)*(errsum)]-denta;
        omegadot = inv(I) * (M - cross(omega, I * omega))+inv(I)     *(disturbance+w);%%% omega is in body frame,
        K_quat = 2; %this enforces the magnitude 1 constraint for the quaternion
        quaterrorl = 1 - (qw^2 + q1^2 + q2^2 + q3^2);
        qdot = -1/2*[0, -p, -q, -r;...
                      p,  0, -r,  q;...
                      q,  r,  0, -p;...
                       r, -q,  p,  0] * quat + K_quat*quaterrorl * quat;
        % Advance system state.
        omega = omega + dt * omegadot;
        quat=quat+qdot*dt; 
        quat=quat/norm(quat);
        quatout(:,ind)=quat;
        errorout(:,ind)=error;
        dentaout(:,ind)=denta;
     
        omegaout(:, ind) =omega;
        inputMout(:, ind) = M;
    end


end
function delta = fC2M( e )
      % delta from the chart definition: Rodrigues Parameters
      delta = [ e; 2 ]/sqrt(4+e'*e);
end    
function G = chartUpdateMatrix( delta )
      G = delta(4)*[ delta(4)   delta(3)  -delta(2);
                    -delta(3)   delta(4)   delta(1);
                     delta(2)  -delta(1)   delta(4)];
end
