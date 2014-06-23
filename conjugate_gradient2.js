compute_E3 = function(A,dim,ro){
    
    var E = new sparse_matrix(A.n,1);
    
    var term1;
    var term2;
    var term3;
    var term4;
    var term5;
    var term6;
    
    for (var i=0; i<dim;i++){
        for (var k=0;k<dim;k++){
            term1 = A.get_element(i*dim+k,i*dim+k);
            
            if (i!=0){
                term2 = A.get_element((i-1)*dim+k, i*dim+k);
                term2 = term2*E.get_element((i-1)*dim+k, (i-1)*dim+k);
                term2 = term2*term2;
            }
            else{
                term2 = 0;
            }
            
            if (k!=0){
                term3 = A.get_element(i*dim+(k-1), i*dim+k);
                term3 = term3*E.get_element(i*dim+(k-1), i*dim+(k-1));
                term3 = term3*term3;
            }
            else{
                term3 = 0;
            }
            if (i!=0 && k!=(dim-1)){
                term4 = A.get_element((i-1)*dim+k, i*dim+k) * A.get_element((i-1)*dim+k, (i-1)*dim+k+1);
                term4 = term4*(E.get_element((i-1)*dim+k, (i-1)*dim+k)*E.get_element((i-1)*dim+k, (i-1)*dim+k));
            }
            else{
                term4 = 0;
            }
            if(k!=0 && i!=(dim-1)){
                term5 = A.get_element(i*dim+(k-1), i*dim+k) * A.get_element(i*dim+(k-1), (i+1)*dim+(k-1));
                term5 = term5*(E.get_element(i*dim+(k-1), i*dim+(k-1)) * E.get_element(i*dim+(k-1), i*dim+(k-1)));
                
            }
            else{
                term5 = 0;
            }
            term6 = term1 -term2 - term3 - ro * (term4+term5);
          

           // document.write("("+i+", "+k+"): "+term1+"; "+term2+"; "+term3+"; "+term4+"; "+term5+"; "+term6+"; ");
            term6 = 1/(Math.sqrt(term6));
            //document.write("("+i+", "+k+"):"+term6+"</br>");
            
            
            E.set_element(i*dim+k, i*dim+k, term6);
            
        }
    }
    
    return E;
}

apply_precon2 = function(A,E,r,dim){
    
    var q = new Array(r.length);
    var z = new Array(r.length);
    
    var term1;
    var term2;
    var term3;
    var term4;
    
    for (var i = 0; i<dim; i++){
    
        for (var k=0; k<dim; k++){
            term1 = r[i*dim+k];
            
            if (i!=0){
                term2 = A.get_element((i-1)*dim+k, i*dim+k) * E.get_element((i-1)*dim+k, (i-1)*dim+k) * q[(i-1)*dim+k];
                //document.write(E.get_element((i-1)*dim+k, (i-1)*dim+k)+"</br>");
            }
            else{
                term2 = 0;
            }
            
            if (k!=0){
                term3 = A.get_element(i*dim+k-1, i*dim+k) * E.get_element(i*dim+k-1, i*dim+k-1) * q[i*dim+k-1];
            }
            else{
                term3 = 0;
            }
            
            term4 = term1 - term2 -term3;
            
            term5 = term4 * E.get_element(i*dim+k, i*dim+k);
            
            q[i*dim+k] = term5;
            //document.write(i+"; "+k+": "+term5+"</br>");
           
            
        }
    }
        for (var i = dim-1; i>-1; i--){
        for (var k=dim-1; k>-1; k--){
            
            term1 = q[i*dim+k];
            
            if (i<(dim-1)){
                term2 = A.get_element(i*dim+k,(i+1)*dim+k) * E.get_element(i*dim+k, i*dim+k) * z[(i+1)*dim+k];
            }
            else{
                term2 = 0;
            }
            
            if (k<(dim-1)){
                term3 = A.get_element(i*dim+k, i*dim+k+1) * E.get_element(i*dim+k, i*dim+k) * z[i*dim+k+1];
            }
            else{
                term3 = 0;
            }
            
            term4 = term1 - term2 - term3;
            term5 = term4 * E.get_element(i*dim+k, i*dim+k);
            z[i*dim+k] = term5;
        }
    }
    q = null;
    return z;
}

function conjugate_gradient2(A,b, iter, dim){
    
    var x = new Array(dim*dim);
    for (var i=0;i<dim*dim;i++){
        x[i]=0;
    }
    var r = b.slice(0);
    var E = compute_E3(A,dim,0.97);
    var roi;
    var roiminus1;
    var p;
    var beta;
    var alpha;
    
    for (var i = 0; i<iter;i++){
        z = apply_precon2(A,E,r,dim);
        roi = numeric.dot(r,z);
        if (i==0){
            p = z.slice(0);
        }
        else{
            beta = roi / roiminus1;
            p = numeric.add(z, numeric.mul(beta,p));
        }
       
        q = A.multiply_by_vector(p);
        
        alpha = roi / (numeric.dot(p,q));
        
        x = numeric.add(x, numeric.mul(alpha,p));

        r = numeric.sub(r, numeric.mul(alpha,q));
        roiminus1 = roi;
        
        check = numeric.norm2(r);
        if (check < 5.0){
            return x;
        }
        
    }
    E = null;
    q = null;
    r = null;
    z = null;
    
    return x;
}