#' Function to define a colormap
#'
#' \code{oColormap} is supposed to define a colormap. It returns a function, which will take an integer argument specifying how many colors interpolate the given colormap.
#'
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{https://html-color-codes.info/color-names/index.html}. It can also be a function of 'grDevices::colorRampPalette'. It can also be "brewer.*" (RColorBrewer palette; see RColorBrewer::display.brewer.all()). It can be colorspace defaults ("rainbow_hcl","heat_hcl","terrain_hcl","diverge_hcl") or other useful ones ("hcl_br","hcl_bp","hcl_bb","hcl_gp","hcl_go","hcl_cp","hcl_cy","hcl_co")
#' @param interpolate use spline or linear interpolation
#' @param data NULL or a numeric vector
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @return 
#' palette.name (a function that takes an integer argument for generating that number of colors interpolating the given sequence) or mapped colors if data is provided.
#' @note The input colormap includes: 
#' \itemize{
#' \item{"jet": jet colormap}
#' \item{"bwr": blue-white-red}
#' \item{"gbr": green-black-red}
#' \item{"wyr": white-yellow-red}
#' \item{"br": black-red}
#' \item{"yr": yellow-red}
#' \item{"wb": white-black}
#' \item{"rainbow": rainbow colormap, that is, red-yellow-green-cyan-blue-magenta}
#' \item{"ggplot2": emulating ggplot2 default color palette}
#' \item{"spectral": emulating RColorBrewer spectral color palette}
#' \item{Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkblue-lightblue-lightyellow-darkorange", "darkgreen-white-darkviolet", "darkgreen-lightgreen-lightpink-darkred". A list of standard color names can be found in \url{https://html-color-codes.info/color-names/index.html}}
#' }
#' @export
#' @seealso \code{\link{oColormap}}
#' @include oColormap.r
#' @examples
#' # 1) define "blue-white-red" colormap
#' palette.name <- oColormap(colormap="bwr")
#' # use the return function "palette.name" to generate 10 colors spanning "bwr"
#' palette.name(10)
#' 
#' # 2) define default colormap from ggplot2
#' palette.name <- oColormap(colormap="ggplot2")
#' # use the return function "palette.name" to generate 3 default colors used by ggplot2
#' palette.name(3)
#' 
#' # 3) define brewer colormap called "RdYlBu"
#' palette.name <- oColormap(colormap="RdYlBu")
#' # use the return function "palette.name" to generate 3 default colors used by ggplot2
#' palette.name(3)
#' 
#' # 4) return mapped colors
#' oColormap(colormap="RdYlBu", data=runif(5))

oColormap <- function(colormap=c("bwr","jet","gbr","wyr","br","yr","rainbow","wb","heat","terrain","topo","cm","ggplot2","jet.top","jet.bottom","jet.both","spectral","ggplot2.top","ggplot2.bottom","ggplot2.both","RdYlBu", "brewer.BrBG","brewer.PiYG","brewer.PRGn","brewer.PuOr","brewer.RdBu","brewer.RdGy","brewer.RdYlBu","brewer.RdYlGn","brewer.Spectral", "brewer.Blues","brewer.BuGn","brewer.BuPu","brewer.GnBu","brewer.Greens","brewer.Greys","brewer.Oranges","brewer.OrRd","brewer.PuBu","brewer.PuBuGn","brewer.PuRd","brewer.Purples","brewer.RdPu","brewer.Reds","brewer.YlGn","brewer.YlGnBu","brewer.YlOrBr","brewer.YlOrRd", "rainbow_hcl","heat_hcl","terrain_hcl","diverge_hcl", "hcl_br","hcl_bp","hcl_bb","hcl_gp","hcl_go","hcl_cp","hcl_cy","hcl_co", "sci_jco","sci_lancet","sci_nejm","sci_locuszoom"), interpolate=c("spline","linear"), data=NULL, zlim=NULL)
{

	interpolate <- match.arg(interpolate)
	
	if(is(colormap,'function')){
		palette.name <- colormap
		
	}else{
	
		if(length(colormap)>1){
			colormap <- colormap[1]
		}
	
		## http://www.cbs.dtu.dk/~eklund/squash/
	
		if(colormap=='ggplot2'){
	
			my_hue_pal <- function (h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1){
				function(n) {
					if ((diff(h)%%360) < 1) {
						h[2] <- h[2] - 360/n
					}
					rotate <- function(x)(x + h.start)%%360 * direction
					hues <- rotate(seq(h[1], h[2], length.out = n))
					grDevices::hcl(hues, c, l)
				}
			}
			palette.name <- my_hue_pal(h=c(0,360)+15, c=100, l=65, h.start=0, direction=1)
		
		}else if(colormap == "jet.top"){
			palette.name <-grDevices::colorRampPalette(c("#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")[-5], interpolate=interpolate)
		}else if(colormap == "jet.bottom"){
			palette.name <-grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F")[-5], interpolate=interpolate)
		}else if(colormap == "jet.both"){
			palette.name <-grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")[c(-1,-9)], interpolate=interpolate)
		
		}else if(colormap == "ggplot2.top"){
			palette.name <-grDevices::colorRampPalette(c("#00C19F","#00B9E3","#619CFF","#DB72FB","#FF61C3"), interpolate=interpolate)
		}else if(colormap == "ggplot2.bottom"){
			palette.name <-grDevices::colorRampPalette(c("#F8766D","#D39200","#93AA00","#00BA38","#00C19F"), interpolate=interpolate)
		}else if(colormap == "ggplot2.both"){
			palette.name <-grDevices::colorRampPalette(c("#F8766D","#D39200","#93AA00","#00BA38","#00C19F","#00B9E3","#619CFF","#DB72FB","#FF61C3")[c(-1,-9)], interpolate=interpolate)
		
		}else if(colormap == "spectral"){
			palette.name <-grDevices::colorRampPalette(rev(c('#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#E6F598','#ABDDA4','#66C2A5','#3288BD')), interpolate=interpolate)

		}else if(colormap == "spectral.both"){
			palette.name <-grDevices::colorRampPalette(rev(c('#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#E6F598','#ABDDA4','#66C2A5','#3288BD'))[c(-1,-9)], interpolate=interpolate)
		}else if(colormap == "spectral.top"){
			palette.name <-grDevices::colorRampPalette(rev(c('#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF')), interpolate=interpolate)
		}else if(colormap == "spectral.bottom"){
			palette.name <-grDevices::colorRampPalette(rev(c('#FFFFBF','#E6F598','#ABDDA4','#66C2A5','#3288BD')), interpolate=interpolate)
		
		}else if(colormap == "RdYlBu"){
			palette.name <-grDevices::colorRampPalette(rev(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695")), interpolate=interpolate)
			
		}else if(grepl("^brewer.",colormap)){
			####################################
			## RColorBrewer::display.brewer.all()
			### Diverging
			#tmp <- rownames(subset(RColorBrewer::brewer.pal.info, category=='div'))
			### Sequential
			#tmp <- rownames(subset(RColorBrewer::brewer.pal.info, category=='seq'))
			#noquote(paste0("brewer.",sort(tmp),collapse='","'))
			####################################
			
			colormap1 <- gsub("^brewer.","",colormap)
			ind <- match(colormap1, rownames(RColorBrewer::brewer.pal.info))
			n <- RColorBrewer::brewer.pal.info$maxcolors[ind]
			if(n==11){
				colors <- rev(RColorBrewer::brewer.pal(n, colormap1))
			}else if(n==9){
				colors <- RColorBrewer::brewer.pal(n, colormap1)
			}
			palette.name <- grDevices::colorRampPalette(colors, interpolate=interpolate)

		}else if(colormap == "rainbow_hcl"){
			#via: noquote(paste0(colorspace::rainbow_hcl(12),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#E495A5","#DB9D85","#C7A76C","#ABB065","#86B875","#5CBD92","#39BEB1","#4CB9CC","#7DB0DD","#ACA4E2","#CD99D8","#E093C3"), interpolate=interpolate)
		}else if(colormap == "heat_hcl"){
			#via: noquote(paste0(colorspace::heat_hcl(12),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#D33F6A","#D95260","#DE6355","#E27449","#E6833D","#E89331","#E9A229","#EAB12A","#E9C037","#E7CE4C","#E4DC68","#E2E6BD"), interpolate=interpolate)
		}else if(colormap == "terrain_hcl"){
			#via: noquote(paste0(colorspace::terrain_hcl(12),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#26A63A","#5FAB25","#84B00E","#A2B309","#BDB723","#D6BA40","#ECBD5D","#FFC07A","#FFC497","#FFC8B3","#FFCFD0","#F1F1F1"), interpolate=interpolate)
		}else if(colormap == "diverge_hcl"){
			#via: noquote(paste0(colorspace::diverge_hcl(12),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#023FA5","#5868AC","#848DBC","#A9AECB","#C8CAD8","#DDDEE0","#E1DDDD","#D9C6C9","#CEA5AC","#BE7E8A","#A94F64","#8E063B"), interpolate=interpolate)
		
		}else if(colormap == "hcl_br"){
			#via: noquote(paste0(colorspace::diverge_hcl(12,c=100,l=c(50,90),power=1),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#4A6FE3","#6D84E1","#8999E1","#A4ADE2","#BDC2E3","#D6D8E3","#E4D4D7","#E6B8C0","#E49CAA","#E18095","#DB627F","#D33F6A"), interpolate=interpolate)
		}else if(colormap == "hcl_bp"){
			#via: noquote(paste0(colorspace::diverge_hcl(12,h=c(255,330),l=c(40,90)),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#155DB5","#617EBE","#8C9CC9","#AEB7D3","#CACEDB","#DEDEE1","#E1DEE0","#DCCAD4","#D3ADC3","#C88AB0","#BA5F9A","#AB1E84"), interpolate=interpolate)
		}else if(colormap == "hcl_bb"){
			#via: noquote(paste0(colorspace::diverge_hcl(12,h=c(246,40),c=96),collapse='","'))
			#palette.name <- grDevices::colorRampPalette(c("#0050B3","#0071B5","#6B92C2","#9DB1CF","#C2CBD9","#DCDEE1","#E0DDDB","#D8C7BE","#CAA995","#B78560","#9E5C00","#7D3200"), interpolate=interpolate)
			palette.name <- grDevices::colorRampPalette(c("#0050B3","#0071B5","#6B92C2","#9DB1CF","#C2CBD9","#DCDEE1","#E0DDDB","#D8C7BE","#CAA995","#B78560","#9E5C00"), interpolate=interpolate)
		}else if(colormap == "hcl_gp"){
			#via: noquote(paste0(colorspace::diverge_hcl(12,h=c(128,330),c=98),l=c(65,90),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#005C0065","#007F0090","#5F9E5F65","#94B99590","#BECFBE65","#DBDFDB90","#E1DDDF65","#DBC4D190","#D1A2BE65","#C377A790","#B33E8E65","#A3007790"), interpolate=interpolate)
		}else if(colormap == "hcl_go"){
			#via: noquote(paste0(colorspace::diverge_hcl(12,h=c(130,43),c=100,l=c(70,90)),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#11C638","#68CF72","#93D699","#B4DBB7","#CDDFCE","#DEE2DF","#E4E0DE","#E8D7CD","#EDCBB4","#F0BC93","#F1AA68","#EF9708"), interpolate=interpolate)
		}else if(colormap == "hcl_cp"){
			#via: noquote(paste0(colorspace::diverge_hcl(12,h=c(180,330),c=59,l=c(75,95)),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#0FCFC0","#75D8CC","#A3DFD7","#C4E6E1","#DCECE9","#EDF0EF","#F1EEF0","#F2E5ED","#F4D7E8","#F6C7E2","#F7B3DC","#F79CD4"), interpolate=interpolate)
		}else if(colormap == "hcl_cy"){
			#via: noquote(paste0(colorspace::diverge_hcl(12,h=c(180,70),c=70,l=c(90,95)),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#22FDEB","#83FAEC","#AEF7EE","#CBF5EF","#E0F3F0","#EDF1F0","#F1F0EE","#F3EEE3","#F6ECD2","#F9E9BD","#FBE5A2","#FEE17F"), interpolate=interpolate)
		}else if(colormap == "hcl_co"){
			#via: noquote(paste0(colorspace::diverge_hcl(12,h=c(180,43),c=70,l=c(90,95)),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#22FDEB","#83FAEC","#AEF7EE","#CBF5EF","#E0F3F0","#EDF1F0","#F2F0EE","#F8ECE5","#FFE8D8","#FFE2C8","#FFDBB3","#FFD49A"), interpolate=interpolate)
		
		}else if(colormap == "sci_jco"){
			#https://cloud.r-project.org/web/packages/ggsci/vignettes/ggsci.html#jco
			#via: noquote(paste0(ggsci::pal_jco("default",alpha=1)(3),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#0073C2FF","#EFC000FF","#868686FF"), interpolate=interpolate)
		}else if(colormap == "sci_nejm"){
			#https://nanx.me/ggsci/reference/pal_nejm.html
			#via: noquote(paste0(ggsci::pal_nejm("default",alpha=1)(8),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF"), interpolate=interpolate)
		}else if(colormap == "sci_lancet"){
			#https://nanx.me/ggsci/reference/pal_lancet.html
			#via: noquote(paste0(ggsci::pal_lancet("lanonc",alpha=1)(6),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"), interpolate=interpolate)
		}else if(colormap == "sci_locuszoom"){
			#https://nanx.me/ggsci/reference/pal_locuszoom.html
			#via: noquote(paste0(ggsci::pal_locuszoom("default",alpha=1)(6),collapse='","'))
			palette.name <- grDevices::colorRampPalette(c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF","#357EBDFF","#9632B8FF"), interpolate=interpolate)
		
		}else if(colormap == "heat"){
			palette.name <- grDevices::heat.colors
		}else if(colormap == "terrain"){
			palette.name <- grDevices::terrain.colors
		}else if(colormap == "topo"){
			palette.name <- grDevices::topo.colors
		}else if(colormap == "cm"){
			palette.name <- grDevices::cm.colors
		
		}else{
			palette.name <- supraHex::visColormap(colormap=colormap)
		}
	}

	########################################
	## return mapped colors
	if(!is.null(data)){	
		if(is.numeric(data)){
			############################
			if(is.null(zlim)){
				vmin <- floor(stats::quantile(data, 0.05))
				vmax <- ceiling(stats::quantile(data, 0.95))
				if(vmin < 0 & vmax > 0){
					vsym <- abs(min(vmin, vmax))
					vmin <- -1*vsym
					vmax <- vsym
				}
				zlim <- c(vmin,vmax)
			}
			data[data<zlim[1]] <- zlim[1]
			data[data>zlim[2]] <- zlim[2]
			############################
			cut_index <- as.numeric(cut(data, breaks=min(data)+(max(data)-min(data))*seq(0, 1, len=64)))
			cut_index[is.na(cut_index)] <- 1
			res <- palette.name(64)[cut_index]
			return(res)
		}
	}
	########################################
		
    invisible(palette.name)
}
