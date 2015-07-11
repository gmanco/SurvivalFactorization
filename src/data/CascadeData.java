package data;


import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

public class CascadeData {
	   
        /*
         * Provides fast access to t_u(c)
         */
        protected SparseDoubleMatrix2D cascadeEvents;
    
        /*
         * For each cascade it contains an ArrayList of cascade events (sorted)
         */
        protected List<CascadeEvent>[] cascadeEventsSrt; 
        
        protected List<WordOccurrence>[] cascadeContent;
         
        /*
         * Set of nodes ids
         */
        protected Set<Integer> nodeSet;
        
        /*
         * Set of cascadesIds
         */
        protected Set<Integer> cascadeSet;
        
        /*
         * Set of word ids
         */
        protected Set<Integer> wordSet;
        
        // properties
        public int n_nodes;
        public int n_cascades;
        public int n_words = 0;
        public double t_max;    
        
        
        public CascadeData(String file_events,String file_content){
            // here we read the file
            // file_cascade is in the format (u,c,t_u(c))
            // file_content is in the format (w,c,n_{w,c})	
            processEventFile(file_events);    
        	
            if (file_content != null)
            		processContentFile(file_content);            
        }
        
        
        
        public List<WordOccurrence> getCascadeContent(int c) {
            return cascadeContent[c];
         }//WordOccurrence
        
        
        /*
         * Return the activation timestamp for the pair
         * (node,cascade). Returns -1 if missing
         */
        public double getActivationTimestamp(int nodeId,int cascadeId){
            Double t=cascadeEvents.get(nodeId, cascadeId);
            if(t==null){
                t=-1.0;
            }
            return t;
        }//getTimestamp
        
       
        
        public List<CascadeEvent> getCascadeEvents(int cascadeId){
            return cascadeEventsSrt[cascadeId];
        }//getSrtEventsForCascade
        
       
        
        
        public void getInfo(){
            System.out.format("Number of nodes: %d\n",n_nodes);
            System.out.format("Number of cascades %d\n",n_cascades);
            System.out.format("Number of words %d\n",n_words);
        }//getInfo        
                
        
        
        public int getNNodes(){
            return n_nodes;
        }
        
        public int getNCascades(){
            return n_cascades;
        }
        
        public int getNWords(){
            return n_words;
        }
                
        public double getTMax(){
            return t_max;
        }
        
       
        public Set<Integer> getNodeIds(){
            return nodeSet;
        }
        
        public Set<Integer> getCascadeIds(){
            return cascadeSet;
        }
        
        public Set<Integer> getWordIds(){
            return wordSet;
        }
        
        /*
         * Read Content file
         */
		private void processContentFile(String file_content) {
			try{
			    wordSet=new TreeSet<Integer>();
			    
			    //read dimensions    
                BufferedReader br=new BufferedReader(new FileReader(file_content));
                String line=br.readLine();
                String tokens[];
                line=br.readLine();//skip header
                
                while(line!=null){
                    tokens=line.split("\t");
                    int word=Integer.parseInt(tokens[0])-1;
                    wordSet.add(word);
                    line=br.readLine();
                }
                br.close();
                
                this.n_words=wordSet.size();
                
                this.cascadeContent=new ArrayList[n_cascades];
               
                br=new BufferedReader(new FileReader(file_content));
                line=br.readLine();
                line=br.readLine();//skip header
             
                List<WordOccurrence>contentTmp;
                while(line!=null){
                    tokens=line.split("\t");
                    int word=Integer.parseInt(tokens[0])-1;
                    int cascadeId=Integer.parseInt(tokens[1])-1;
                    int cnt=Integer.parseInt(tokens[2]);
                   
                    contentTmp=cascadeContent[cascadeId];
                    if(contentTmp==null){
                        contentTmp=new ArrayList<WordOccurrence>();
                        cascadeContent[cascadeId]=contentTmp;
                    }
                    contentTmp.add(new WordOccurrence(word, cnt));
                    
                   
                    line=br.readLine();
                }
                br.close();
			    
			}catch(Exception e){
			    e.printStackTrace();
			}
			
		}//processContentFile


		/*
		 * Read event file
		 */
		private void processEventFile(String file_events) {
		    try{
			
		        this.nodeSet=new TreeSet<Integer>();
	            this.cascadeSet=new TreeSet<Integer>();
	            this.t_max=-1;
	            
    		    //read dimensions    
    		    BufferedReader br=new BufferedReader(new FileReader(file_events));
    			String line=br.readLine();
    			String tokens[];
    			line=br.readLine();//skip header
    			
    			while(line!=null){
    			    tokens=line.split("\t");
    			    int nodeId=Integer.parseInt(tokens[0])-1;
    			    int cascadeId=Integer.parseInt(tokens[1])-1;
    			    nodeSet.add(nodeId);
    			    cascadeSet.add(cascadeId);
    			    line=br.readLine();
    			}
    			br.close();
			
    			//init variables
    			this.n_nodes=nodeSet.size();
    			this.n_cascades=cascadeSet.size();
    			this.cascadeEvents=new SparseDoubleMatrix2D(n_nodes,n_cascades);
    			this.cascadeEventsSrt=new ArrayList[n_cascades];
    			
    			//read cascades
    			 br=new BufferedReader(new FileReader(file_events));
    	         line=br.readLine();
    	         line=br.readLine();//skip header
    	        
    	         List<CascadeEvent>tmp;
    	         while(line!=null){
    	                tokens=line.split("\t");
    	                int nodeId=Integer.parseInt(tokens[0]) -1 ;
    	                int cascadeId=Integer.parseInt(tokens[1]) -1;
    	                double timestamp=Double.parseDouble(tokens[2]);
    	              
    	                //set the timestamp
    	                cascadeEvents.set(nodeId, cascadeId, timestamp);
    	                
    	                // add event to CascadeEventsSrt
    	                tmp=cascadeEventsSrt[cascadeId];
    	                if(tmp==null){
    	                    tmp=new ArrayList<CascadeEvent>();
    	                    cascadeEventsSrt[cascadeId]=tmp;
    	                }
    	                tmp.add(new CascadeEvent(nodeId,timestamp));
    	              
    	                //check t_max
    	                if(this.t_max<timestamp){
    	                    this.t_max=timestamp;
    	                }
    	                
    	          line=br.readLine();     
    	         }
    	         
    	         br.close();
    	         
    	         //sort events within cascades
    	         for(int cascadeId:cascadeSet){
                     tmp=cascadeEventsSrt[cascadeId];
                     Collections.sort(tmp);        
                 }//sort 
    	         
    		    }catch(Exception e ){
    		        e.printStackTrace();
    		    }
		}//processEventFile

      public static void main(String[] args) throws Exception{
        String file_events="/Users/barbieri/Dropbox/shared ICAR/SurvivalFactorization/exp/meme_tracker/debug/cleaned_debug_activations";
        String file_content="/Users/barbieri/Dropbox/shared ICAR/SurvivalFactorization/exp/meme_tracker/debug/cleaned_debug_content";
    
        CascadeData data=new CascadeData(file_events, file_content);
        data.getInfo();
      
      }



   
        

}
