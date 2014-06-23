Array.prototype.insert = function (index, item){
    this.splice(index, 0, item);
};

function sparse_matrix(n, expected_nonzeros_per_row){
    
    this.n = n;
    this.index = new Array(n);
    this.value = new Array(n);
    this.number_of_values_in_row = new Array(n);
    
    for (var i = 0;i<n;i++){
        this.number_of_values_in_row[i]=0;
    }
    
    for (var i = 0;i<n;i++){
        this.index[i] = new Array(expected_nonzeros_per_row);
        this.value[i] = new Array(expected_nonzeros_per_row);
    }
    
    this.set_element = function (i,j,value){
        
        for (var k=0;k<this.number_of_values_in_row[i];k++){
        
            if (this.index[i][k]==j){
                this.value[i][k]=value;
                return;
            }
            else if(this.index[i][k]>j){
                this.index[i].insert(k,j);
                this.value[i].insert(k,value);
                this.number_of_values_in_row[i]++;
                return;
            }
        }
        
        this.index[i][this.number_of_values_in_row[i]]=j;
        this.value[i][this.number_of_values_in_row[i]]=value;
        this.number_of_values_in_row[i]++;
        
        
    }
    
    this.get_element = function(i,j){
        
        for (var k=0; k<this.number_of_values_in_row[i]; k++){
            if (this.index[i][k]==j){
                return this.value[i][k];
            }
            else if (this.index[i][k]>j){
                return 0;
            }
        }
        return 0;
        
    }
    
    this.print_matrix = function(){

        
        var counter;
        for (var i = 0; i<this.index.length; i++){
            counter=0;
            for (var k = 0; k<this.index.length;k++){
                
                if (counter<this.number_of_values_in_row[i]){
               
                    if (this.index[i][counter]==k){
                        document.write(this.value[i][counter]);
                        counter++;
                    }
                    else{
                        document.write("0");
                    }
                }
                else{
                    document.write("0");
                }
                
                if (k==this.index.length-1){
                    document.write("</br>");
                }
                else{
                    document.write(", ");
                }
            }
        }
        
    }
    
    this.multiply_by_vector = function(b){
            
            var result = new Array(b.length);
            for (var i=0;i<b.length;i++){
                result[i]=0;
                for (var k=0;k<this.number_of_values_in_row[i]; k++){
                    result[i]+=this.value[i][k]*b[this.index[i][k]];
                }
            }
            
            return result;
            
        }
        
    
}