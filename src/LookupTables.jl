# Copyright (C) 2016 Bernat Zaragoza Travieso
# See accompanying license file.

module LookupTables


export LookupTable2d
export lookup


# Rows vector: vertical vector containing the values for each row
# Columns vector: horizontal vector containing the values for each colum

struct LookupTable2d{T ,U ,V }
	rowsvector::Array{T,1}
	columnsvector::Array{U,1}
	output::Array{V,2}
	
	function LookupTable2d(rowsvector, columnsvector, output)
		# arguments dimensions check
		inputcolumns = length(columnsvector)
		inputrows = length(rowsvector)
		outputcolumns = size(output,2)
		outputrows = size(output,1)
		rowerror, columnerror = false, false
		if inputrows != outputrows && inputcolumns != outputcolumns
			error("Invalid arguments for LookupTable2d constructor: number of rows and number of columns don't match")
		elseif inputcolumns != outputcolumns
			error("Invalid arguments for LookupTable2d constructor: number of columns don't match")
		elseif inputrows != outputrows
			error("Invalid arguments for LookupTable2d constructor: number of rows don't match")
		elseif inputcolumns==inputrows==1
			error("Invalid arguments for LookupTable2d constructor: both table dimensions must be > 1")
		end
		
		# check that rows vector and columns vector are strictly monotonic
		for i =2:inputrows
			if rowsvector[i] <= rowsvector[i-1]
				rowerror = true
				break
			end
		end
		for i =2:inputcolumns
			if columnsvector[i] <= columnsvector[i-1]
				columnerror = true
				break
			end
		end
		if rowerror && columnerror
			error("Provided rows vector and columns vector must be strictly monotonic")
		elseif rowerror
			error("Provided rows vector must be strictly monotonic")
		elseif columnerror
			error("Provided columns vector must be strictly monotonic")
		end
		new{typeof(rowsvector[1]), typeof(columnsvector[1]), typeof(output[1])}(rowsvector,columnsvector,output)
	end
end


#= Returns the index of the closest smallest number to "value" from the array "l". 
   If "value" exists in "l", it returns the index of this value.
   If "value" is out of bounds, it returns 1 (for a smaller "value")
   or the length-1 of "l" (for a higher "value") 
   Ex.		l= [1,4,5,7,8]
   rangesearch(l, 4) -> low = 2
   rangesearch(l, 6) -> low = 3
   rangesearch(l, 10) -> low = 4
   
 =#
function rangesearch(l, value)
	low = 1
	high = length(l)
	while low+1 < high
			mid = round(Int,(low+high)/2)
			if l[mid] > value 
				high = mid
			elseif l[mid] <= value 
				low = mid
			end
	end		 
	return low
end


#= method:	0 Interpolation-Extrapolation
			1 Interpolation-Use End Values
			2 Use Input Nearest
			3 Use Input Below
			4 Use Input Above
=#

function lookup(row::Real, column::Real, table::LookupTable2d, method::Integer=0)
	# Binary search in rows vector
	rowlow = rangesearch(table.rowsvector, row)
	columnlow = rangesearch(table.columnsvector, column)
	# Behaviour according to method
	if method == 0
		rowhigh = rowlow + 1
		columnhigh = columnlow + 1
		# Calculate column and row factors
		rowfactor = (row-table.rowsvector[rowlow])/(table.rowsvector[rowhigh]-table.rowsvector[rowlow])
		columnfactor = (column-table.columnsvector[columnlow])/(table.columnsvector[columnhigh]-table.columnsvector[columnlow])	
		# For row low
		outputcolumnΔlow = table.output[rowlow,columnhigh]-table.output[rowlow,columnlow]
		# For row high
		outputcolumnΔhigh = table.output[rowhigh,columnhigh]-table.output[rowhigh,columnlow]
		# For this column
		outputrowΔ = table.output[rowhigh,columnlow] + outputcolumnΔhigh*columnfactor - (table.output[rowlow,columnlow] + outputcolumnΔlow*columnfactor)
		return table.output[rowlow,columnlow] + columnfactor*outputcolumnΔlow + rowfactor*outputrowΔ	
	elseif method == 1
		rowhigh = rowlow + 1
		columnhigh = columnlow + 1		
		if column >= table.columnsvector[columnhigh]
			columnfactor =  0
			outputcolumnΔlow = 0
			columnlow = columnhigh
		elseif column <= table.columnsvector[columnlow]
			columnfactor =  0
			outputcolumnΔlow = 0
		else
			# For row low
			outputcolumnΔlow = table.output[rowlow,columnhigh]-table.output[rowlow,columnlow]
		end
		if row >= table.rowsvector[rowhigh]
			rowfactor =  0
			outputrowΔ = 0
			rowlow = rowhigh
		elseif row <= table.rowsvector[rowlow]
			rowfactor =  0
			outputrowΔ = 0
		else
			# For row high
			outputcolumnΔhigh = table.output[rowhigh,columnhigh]-table.output[rowhigh,columnlow]
			# Calculate row factor
			rowfactor = (row-table.rowsvector[rowlow])/(table.rowsvector[rowhigh]-table.rowsvector[rowlow])
			# For this column
			outputrowΔ = table.output[rowhigh,columnlow] + outputcolumnΔhigh*columnfactor - (table.output[rowlow,columnlow] + outputcolumnΔlow*columnfactor)
		end
		return table.output[rowlow,columnlow] + columnfactor*outputcolumnΔlow + rowfactor*outputrowΔ
	elseif method == 2
		rowhigh = rowlow + 1
		columnhigh = columnlow + 1
		if row < (table.rowsvector[rowlow]+table.rowsvector[rowhigh])/2
			nearestrow = rowlow
		else
			nearestrow = rowhigh
		end
		if column < (table.columnsvector[columnlow]+table.columnsvector[columnhigh])/2
			nearestcolumn = columnlow
		else
			nearestcolumn = columnhigh
		end
		return table.output[nearestrow,nearestcolumn]
	elseif method == 3
		return table.output[rowlow,columnlow]
	elseif method == 4
		if table.rowsvector[rowlow] == table.rowsvector[row]
			rowhigh = rowlow
		else
			rowhigh = rowlow + 1
		end
		if table.columnsvector[columnlow] == table.columnsvector[column]
			columnhigh = columnlow
		else
			columnhigh = columnlow + 1
		end
		return table.output[rowhigh, columnhigh]
	else
		error("Invalid method. Use: (0) Interpolation-Extrapolation , (1) Interpolation-Use End Values, (2) Use Input Nearest, (3) Use Input Below or (4) Use Input Above")
	end
end


end #module
