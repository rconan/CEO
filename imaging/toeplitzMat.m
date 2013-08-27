classdef toeplitzMat
    %% TOEPLITZMAT Defines a toeplitz matrix object
    % 
    % T = toeplitzMat(col,row) creates a toeplitz matrix object with col as
    % the leftmost column and row as the top row
    % 
    % T = toeplitzMat(n,m,elements) creates a nxm toeplitz matrix object
    % with the leftmost column givem by elements(m+1:end) and the top row
    % given by elements(1:m)
   
    properties
       nRow;
       nCol;
       elements;
       bytes;
    end
    
    methods
        
        function obj = toeplitzMat(varargin)
            switch nargin
                case 1
                    T = varargin{1};
                    obj.nRow = length(T.elements);
                    obj.nCol = length(T.elements);
                    a = arrayfun( @(x) transpose(x) , T.elements(end-1:-1:1) , 'uniformOutput', false);
                    obj.elements = [T.elements;cat(1,a{:})];
                case 2
                    col = varargin{1};
                    col = col(:);
                    row = varargin{2};
                    row = row(:);
                    if isnumeric(col) && (col(1)~=row(1))
                        error('toeplitzMat:toeplitzMat','First number of row and column must match!')
                    end
                    obj.nRow = length(col);
                    obj.nCol = length(row);
                    obj.elements = [ row(obj.nCol:-1:2) ; col ];
                case 3
                    obj.nRow = varargin{1};
                    obj.nCol = varargin{2};
                    obj.elements = varargin{3};
                    if strcmp(obj,'toeplitzMat') && (length(obj.elements)~=obj.nRow+obj.nCol-1)
                        error('toeplitzMat:toeplitzMat','Elements must contain %d values!',obj.nRow+obj.nCol-1)
                    end
                    obj.elements = obj.elements(:);
                otherwise                    
                    help('toeplitzMat')
                    error('toeplitzMat:toeplitzMat','Wrong constructor call!')
            end
            if isnumeric(obj.elements)
                obj.bytes = length(obj.elements)*8;
            else
                obj.bytes = sum( [ obj.elements.bytes ] );
            end
        end
        
%         function display(obj)
%             nObj = length(obj);
%             if nObj>1
%                 fprintf('%d(%dx%d) Toeplitz matrix: %d elements stored\n',...
%                     nObj,obj(1).nRow,obj(1).nCol,length(obj(1).elements));
%             else
%                 fprintf('%dx%d Toeplitz matrix: %d elements stored\n',...
%                     obj.nRow,obj.nCol,length(obj.elements));
%             end
%         end
        
        function t = full(obj)
            %% FULL Returns the full matrix
            %
            % FT = full(T) return a full Toeplitz matrix from the
            % toeplitzMat object T
            
            p = obj.nCol;
            m = obj.nRow;
            x = obj.elements;                 % build vector of user data
            cidx = (0:m-1)';
            ridx = p:-1:1;
            t = cidx(:,ones(p,1)) + ridx(ones(m,1),:);  % Toeplitz subscripts
            t = x(t);                                   % actual data
            if isa(t(1),'toeplitzMat')
                 t = cell2mat(...
                     arrayfun( @(x) full(x) , t,'uniformOutput',false));
            end
        end
        
        function obj = transpose(obj)
            %% TRANSPOSE Toeplitz matrix transpose
            %
            % TT = transpose(T) returns the transpose of the Toeplitz
            % matrix into a new toeplitzMat object
            
            obj.elements = flipud(obj.elements);
            if isa(obj.elements(1),'toeplitzMat')
                a = arrayfun( @(x) transpose(x) , obj.elements , 'uniformOutput', false);
                obj.elements = cat(1 , a{:} );
            end
        end
        
        function obj = plus(obj1,obj2)
            if (obj1.nRow==obj2.nRow) && (obj1.nCol==obj2.nCol)
                obj = toeplitzMat(obj1.nRow,obj1.nCol,obj1.elements+obj2.elements);
            else
                error('toeplitzMat:plus','Both Toeplitz matrix must be the same size!')
            end
        end
            
    end   
    
    methods (Static)
        
        function demo
            
            n = 4;
            T = toeplitzMat(n,n,1:2*n-1);
            disp(T)
            full(T)
            
            m = 7;
            T = toeplitzMat(n,m,1:n+m-1);
            disp(T)
            full(T)
            
            nB = 5;
            for k=1:2*nB-1
                T(k) = toeplitzMat(n,n,(1:2*n-1)+rand*n*2);
            end
            TBT = toeplitzMat(nB,nB,T);
            disp(TBT)
            FTBT = full(TBT);
            figure,imagesc(FTBT),axis equal tight
            
            clear T
            for k=1:2*n-1
                T(k) = toeplitzMat(n,m,(1:n+m-1)+rand*n);
            end
            TBT = toeplitzMat(n,n,T);
            disp(TBT)
            FTBT = full(TBT);
            figure,imagesc(FTBT),axis equal tight
            
            nBB = 3;
            nB = 5;
            for kB=1:2*nBB-1
                for k=1:2*nB-1
                    T(k) = toeplitzMat(n,n,(1:2*n-1)+rand*n*2);
                end
                TBT(kB) = toeplitzMat(nB,nB,T);
            end
            TBTBT = toeplitzMat(nBB,nBB,TBT);
            disp(TBTBT)
            FTBTBT = full(TBTBT);
            figure,imagesc(FTBTBT),axis equal tight
        end
        
    end
    
end