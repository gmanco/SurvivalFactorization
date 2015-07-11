package it.cnr.adalab.surivalfactorization;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

public class CascadeData {
	   
        double[][][] CascadeEvents; // is a set in which each element is a Cascade (2-dim matrix) user x time
        SparseDoubleMatrix2D CascadeContent; // is a 2-dim sparse matrix word x cascade
        
        SparseDoubleMatrix2D UserCascade_Matrix;
        
        int n_users;
        int n_cascades;
        int n_words = 0;
        double t_max;    
        
        public CascadeData(String file_events,String file_content){
            // here we read the file
            // file_cascade is in the format (u,c,t_u(c))
            // file_content is in the format (w,c,n_{w,c})
        		processEventFile(file_events);    
        	
            if (file_content != null)
            		processContentFile(file_content);            
        }
        


		private void processContentFile(String file_content) {
			// TODO Auto-generated method stub
			
		}



		private void processEventFile(String file_events) {
			// TODO Auto-generated method stub
			
		}



        public void getInfo(){
			System.out.format("Number of users: %d\n",n_users);
			System.out.format("Number of cascades %d\n",n_cascades);
					System.out.format("Number of words %d\n",n_words);
		}        
        
               
        public double getN_Users(){
            return n_users;
        }
        
        public double getN_Cascades(){
            return n_cascades;
        }
        
        public double getN_Words(){
            return n_words;
        }
                
        public double getT_Max(){
            return t_max;
        }
        
        public SparseDoubleMatrix2D getCascadeContent(){
            return CascadeContent;
        }
        
        public double[][] getCascadeEvents(int c){
                return CascadeEvents[c];
        }
        
        public HashMap<Integer,Double> getContent(int c){
        		return CascadeContent.getColumn(c);
        }        
        
        public Set<Integer> getUsersForCascade(int c){
                return CascadeContent.getColumn(c).keySet();
        }        
        

}
