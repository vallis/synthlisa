<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns="http://www.w3.org/1999/xhtml" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="XSIL">
    <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
        <head>
            <link href="lisa-xml.css" rel="Stylesheet" type="text/css" media="screen"/>    
            <title>Mock LISA Data Challenge XML File Format, v. 1.0</title>
        </head>
        <body>
            <div class="FileInfo">
                <h2>Mock LISA Data Challenge XML File Format, v. 1.0</h2>

                <div class="FileData">   
                    <h3>File Info</h3>

                    <table>
        
                    <tr>
                        <td><b>Authors</b></td>
                        <td><xsl:for-each select="Param[@Name='Author']">
                            <xsl:value-of select="."/><br/>
                        </xsl:for-each></td>
                        <td></td>
                    </tr>
        
                    <xsl:for-each select="Param[@Name!='Author']">
                        <tr>
                            <td><b><xsl:value-of select="./@Name"/></b></td>
                            <td><xsl:value-of select="."/></td>
                            <td>
                                <xsl:if test="@Unit != ''">
                                    <xsl:value-of select="./@Unit"/>
                                </xsl:if>
                                <xsl:if test="@Type != ''">
                                    <xsl:value-of select="./@Type"/>
                                </xsl:if>
                            </td>
                        </tr>
                    </xsl:for-each>
        
                    </table>
                    
                    <xsl:for-each select="Comment">
                        <xsl:call-template name="handle-comment"/>
                    </xsl:for-each>
                </div>
            </div>

            <xsl:for-each select="XSIL[@Type='LISACode']">
                <div class="LISAData">                
                    <h2>LISACode parameters</h2>

					<div class="FileData">
						<xsl:call-template name="handle-parameter-table"/>

						<xsl:for-each select="XSIL">
							<h3><xsl:value-of select="./@Name"/> (<xsl:value-of select="./@Type"/>)</h3>

							<div class="TimeSeries">
								<xsl:call-template name="handle-parameter-table"/>
							</div>
						</xsl:for-each>

	                    <xsl:for-each select="Comment">
	                        <xsl:call-template name="handle-comment"/>
	                    </xsl:for-each>
					</div>
                </div>
            </xsl:for-each>

            <xsl:for-each select="XSIL[@Type='LISAData']">
                <div class="LISAData">                
                    <h2>LISA data</h2>
                    
                    <xsl:for-each select="XSIL">
                        <div class="LISA">
                            <h3><xsl:value-of select="./@Name"/> (<xsl:value-of select="./@Type"/>)</h3>

                            <xsl:call-template name="handle-parameter-table"/>

                            <xsl:for-each select="XSIL[@Type='TimeSeries']">
                                <xsl:call-template name="handle-timeseries"/>
                            </xsl:for-each>                           

                            <xsl:for-each select="Comment">
                                <xsl:call-template name="handle-comment"/>
                            </xsl:for-each>
                        </div>
                    </xsl:for-each>
                    
                    <xsl:for-each select="Comment">
                        <xsl:call-template name="handle-comment"/>
                    </xsl:for-each>
                </div>
            </xsl:for-each>
            
            <xsl:for-each select="XSIL[@Type='SourceData']">
                <div class="SourceData">                
                    <h2>Source data</h2>
                    
                    <xsl:for-each select="XSIL[@Type='SampledPlaneWave' or @Type='PlaneWave']">
                        <div class="Wave">
                            <h3><xsl:value-of select="./@Name"/> (<xsl:value-of select="./@Type"/>)</h3>

                            <xsl:call-template name="handle-parameter-table"/>

                            <xsl:for-each select="XSIL[@Type='TimeSeries']">
                                <xsl:call-template name="handle-timeseries"/>
                            </xsl:for-each>                           

                            <xsl:for-each select="Comment">
                                <xsl:call-template name="handle-comment"/>
                            </xsl:for-each>
                        </div>
                    </xsl:for-each>

                    <xsl:for-each select="XSIL[@Type='PlaneWaveTable']">
                        <div class="Wave">
                            <h3><xsl:value-of select="./@Name"/> (<xsl:value-of select="./@Type"/>)</h3>

                            <xsl:call-template name="handle-parameter-table"/>

                            <xsl:for-each select="Table">
                                <xsl:call-template name="handle-table-columns"/>
                            </xsl:for-each>

                            <xsl:for-each select="XSIL[@Type='TimeSeries']">
                                <xsl:call-template name="handle-timeseries"/>
                            </xsl:for-each>

                            <xsl:for-each select="Comment">
                                <xsl:call-template name="handle-comment"/>
                            </xsl:for-each>
                        </div>
                    </xsl:for-each>
                    
                    <xsl:for-each select="Comment">
                        <xsl:call-template name="handle-comment"/>
                    </xsl:for-each>
                </div>
            </xsl:for-each>
                   
            <xsl:for-each select="XSIL[@Type='TDIData']">
                <div class="TDIData">                
                    <h2>TDI data</h2>
                    
                    <xsl:for-each select="XSIL">
                        <div class="Observable">
                            <h3><xsl:value-of select="./@Name"/> (<xsl:value-of select="./@Type"/>)</h3>

                            <xsl:call-template name="handle-parameter-table"/>

                            <xsl:for-each select="XSIL[@Type='TimeSeries']">
                                <xsl:call-template name="handle-timeseries"/>
                            </xsl:for-each>

                            <xsl:for-each select="Comment">
                                <xsl:call-template name="handle-comment"/>
                            </xsl:for-each>
                        </div>
                    </xsl:for-each>
                    
                    <xsl:for-each select="Comment">
                        <xsl:call-template name="handle-comment"/>
                    </xsl:for-each>
                </div>
            </xsl:for-each>
            
            <xsl:for-each select="XSIL[@Type='NoiseData']">
                <div class="NoiseData">                
                    <h2>LISA noise data</h2>
                    
                    <xsl:for-each select="XSIL">
                        <div class="Noise">
                            <h3><xsl:value-of select="./@Name"/> (<xsl:value-of select="./@Type"/>)</h3>

                            <xsl:call-template name="handle-parameter-table"/>

                            <xsl:for-each select="XSIL[@Type='TimeSeries']">
                                <xsl:call-template name="handle-timeseries"/>
                            </xsl:for-each>                           

                            <xsl:for-each select="Comment">
                                <xsl:call-template name="handle-comment"/>
                            </xsl:for-each>
                        </div>
                    </xsl:for-each>

                    <xsl:for-each select="XSIL[@Type='PseudoRandomNoiseTable']">
                        <div class="Noise">
                            <h3><xsl:value-of select="./@Name"/> (<xsl:value-of select="./@Type"/>)</h3>

                            <xsl:call-template name="handle-parameter-table"/>

                            <xsl:for-each select="Table">
                                <xsl:call-template name="handle-table-columns"/>
                            </xsl:for-each>

                            <xsl:for-each select="XSIL[@Type='TimeSeries']">
                                <xsl:call-template name="handle-timeseries"/>
                            </xsl:for-each>

                            <xsl:for-each select="Comment">
                                <xsl:call-template name="handle-comment"/>
                            </xsl:for-each>
                        </div>
                    </xsl:for-each>
                    
                    <xsl:for-each select="Comment">
                        <xsl:call-template name="handle-comment"/>
                    </xsl:for-each>
                </div>
            </xsl:for-each>
            
        </body>
    </html>
</xsl:template>

<xsl:template name="handle-parameter-table">
    <table>                            
        <xsl:for-each select="Param">
            <tr>
                <td><xsl:value-of select="./@Name"/></td>
                <td><xsl:value-of select="."/></td>
                <td><xsl:if test="@Unit != ''">
                    <xsl:value-of select="./@Unit"/>
                </xsl:if></td>
            </tr>
        </xsl:for-each>
    </table>
</xsl:template>

<xsl:template name="handle-table-columns">
    <table>
        <tr>
            <td><b>Parameter Table</b></td>
            <xsl:if test="Stream/@Type = 'Remote'">
                <td>Filename</td>
                <td><a href="{Stream}"><xsl:value-of select="Stream"/></a></td>                     
            </xsl:if>
            <xsl:if test="Stream/@Type = 'Local'">
                <td>(see below)</td><td></td>
            </xsl:if>
        </tr>
        <tr>
            <td>(<xsl:value-of select="Dim[@Name='Length']"/> rows x <xsl:value-of select="Dim[@Name='Records']"/> columns)</td>
            <td>Encoding</td>
            <td><xsl:value-of select="Stream/@Encoding"/></td>
        </tr>
        <tr>
            <td></td>
            <td>Type</td>
            <td><xsl:value-of select="Stream/@Type"/></td>
        </tr>        
        <xsl:for-each select="Column">
            <tr>
                <td>column <xsl:value-of select="position()"/>:</td>
                <td><xsl:value-of select="./@Name"/> (<xsl:value-of select="./@Type"/>)</td>
                <td><xsl:if test="@Unit != ''">
                    <xsl:value-of select="./@Unit"/>
                </xsl:if></td>
            </tr>
        </xsl:for-each>
    </table>
	
	<xsl:if test="Stream/@Type = 'Local'">
		<pre>
        	<xsl:call-template name="pre-text">
        		<xsl:with-param name="string" select="Stream" />
        	</xsl:call-template>
		</pre>
	</xsl:if>
</xsl:template>

<xsl:template name="pre-text">
   <xsl:param name="string" select="." />
   <xsl:choose>
      <xsl:when test="contains($string, '&#xA;')">
         <xsl:value-of select="normalize-space(substring-before($string, '&#xA;'))" /><xsl:text>
</xsl:text><xsl:call-template name="pre-text">
            <xsl:with-param name="string" select="substring-after($string, '&#xA;')" />
         </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
         <xsl:value-of select="normalize-space($string)" />
      </xsl:otherwise>
   </xsl:choose>
</xsl:template>

<xsl:template name="add-line-breaks">
   <xsl:param name="string" select="." />
   <xsl:choose>
      <xsl:when test="contains($string, '&#xA;')">
         <xsl:value-of select="substring-before($string, '&#xA;')" />
         <br/>
         <xsl:call-template name="add-line-breaks">
            <xsl:with-param name="string" select="substring-after($string, '&#xA;')" />
         </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
         <xsl:value-of select="$string" />
      </xsl:otherwise>
   </xsl:choose>
</xsl:template>

<xsl:template name="handle-timeseries">
    <div class="TimeSeries">
        <h4>TimeSeries: <xsl:value-of select="./@Name"/></h4>

        <table>
            <xsl:for-each select="Param">
                <tr>
                    <td><xsl:value-of select="./@Name"/></td>
                    <td><xsl:value-of select="."/></td>
                    <td><xsl:if test="@Unit != ''">
                        <xsl:value-of select="./@Unit"/>
                    </xsl:if></td>
                </tr>
            </xsl:for-each>
            <xsl:for-each select="Array">
                    <tr>
                        <td><b>Array Stream: <xsl:value-of select="./@Name"/></b></td>
                        <td>Filename</td>
                        <td><a href="{Stream}"><xsl:value-of select="./Stream[@Type='Remote']"/></a></td>
                    </tr><tr>
                        <td></td>
                        <td>Encoding</td>                                         
                        <td><xsl:value-of select="./Stream[@Type='Remote']/@Encoding"/></td>
                    </tr><tr>
                        <td></td>
                        <td>Type</td>                                         
                        <td><xsl:value-of select="./@Type"/></td>
                    </tr><tr>
                        <td></td>
                        <td>Unit</td>                                         
                        <td><xsl:value-of select="./@Unit"/></td>                                            
                    </tr><tr>
                        <td></td>
                        <td>Length</td>
                        <td><xsl:value-of select="./Dim[@Name='Length']"/></td>
                    </tr><tr>
                        <td></td>
                        <td>Records</td>
                        <td><xsl:value-of select="./Dim[@Name='Records']"/></td>                                            
                    </tr>
            </xsl:for-each>
        </table>
        
        <xsl:for-each select="Comment">
            <xsl:call-template name="handle-comment"/>
        </xsl:for-each>
    </div>
</xsl:template>

<xsl:template name="handle-comment">
    <div class="Comment"><p>
        <xsl:call-template name="add-line-breaks">
            <xsl:with-param name="string">
                <xsl:value-of select="."/>
            </xsl:with-param>
        </xsl:call-template>
    </p></div>
</xsl:template>
                        
</xsl:stylesheet>
